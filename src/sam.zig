const std = @import("std");
const testing = std.testing;
const Allocator = std.mem.Allocator;

// Read-only SAM text parser.
// Reference: https://samtools.github.io/hts-specs/SAMv1.pdf
//
// Usage:
//   var r: std.Io.Reader = .fixed(bytes);
//   var sr = sam.reader(allocator, &r);
//   var header = try sr.readHeader();
//   defer header.deinit(allocator);
//   while (try sr.next(&header)) |record| { ... }
//
// The `Header` owns its strings and is valid until `deinit`. An `Alignment`
// returned by `next` holds slices into the reader's line buffer and is only
// valid until the following `next` call.

pub const SortOrder = enum { unknown, unsorted, query_name, coordinate };

pub const GroupAlignment = enum {
    none, // default
    query, // alignments are grouped by QNAME
    reference, // alignments are grouped by RNAME/POS
};

pub const SamAlignmentFlags = struct {
    const Self = @This();

    paired: bool = false, // 0x1 the read is paired in sequencing
    proper_pair: bool = false, // 0x2 each segment properly aligned according to the aligner
    unmap: bool = false, // 0x4 segment is unmapped
    munmap: bool = false, // 0x8 next segment in the template unmapped
    reverse: bool = false, // 0x10 SEQ being reverse complemented
    mate_reverse: bool = false, // 0x20 SEQ of the next segment being reverse complemented
    read1: bool = false, // 0x40 first segment in the template
    read2: bool = false, // 0x80 last segment in the template
    secondary: bool = false, // 0x100 secondary alignment
    qc_failure: bool = false, // 0x200 not passing filters (QC)
    dup: bool = false, // 0x400 PCR or optical duplicate
    sup_alignment: bool = false, // 0x800 supplementary alignment

    pub const PairedBit = 0x1;
    pub const ProperPairedBit = 0x2;
    pub const Unmapped = 0x4;
    pub const MateUnmapped = 0x8;
    pub const Reverse = 0x10;
    pub const MateReverse = 0x20;
    pub const Read1 = 0x40;
    pub const Read2 = 0x80;
    pub const Secondary = 0x100;
    pub const QCFailure = 0x200;
    pub const Dup = 0x400;
    pub const SupplementaryAlignment = 0x800;

    pub fn toSamAlignmentFlags(value: u16) SamAlignmentFlags {
        return .{
            .paired = (PairedBit & value) != 0,
            .proper_pair = (ProperPairedBit & value) != 0,
            .unmap = (Unmapped & value) != 0,
            .munmap = (MateUnmapped & value) != 0,
            .reverse = (Reverse & value) != 0,
            .mate_reverse = (MateReverse & value) != 0,
            .read1 = (Read1 & value) != 0,
            .read2 = (Read2 & value) != 0,
            .secondary = (Secondary & value) != 0,
            .qc_failure = (QCFailure & value) != 0,
            .dup = (Dup & value) != 0,
            .sup_alignment = (SupplementaryAlignment & value) != 0,
        };
    }

    pub fn toValue(self: Self) u16 {
        var v: u16 = 0;
        if (self.paired) v |= PairedBit;
        if (self.proper_pair) v |= ProperPairedBit;
        if (self.unmap) v |= Unmapped;
        if (self.munmap) v |= MateUnmapped;
        if (self.reverse) v |= Reverse;
        if (self.mate_reverse) v |= MateReverse;
        if (self.read1) v |= Read1;
        if (self.read2) v |= Read2;
        if (self.secondary) v |= Secondary;
        if (self.qc_failure) v |= QCFailure;
        if (self.dup) v |= Dup;
        if (self.sup_alignment) v |= SupplementaryAlignment;
        return v;
    }
};

// CIGAR ---------------------------------------------------------------------

/// CIGAR operations, ordered to match the "MIDNSHP=XB" character table.
pub const CigarOp = enum(u8) {
    match, // M
    ins, // I
    del, // D
    ref_skip, // N
    soft_clip, // S
    hard_clip, // H
    pad, // P
    equal, // =
    diff, // X
    back, // B

    pub const chars = "MIDNSHP=XB";

    pub fn fromChar(c: u8) ?CigarOp {
        return switch (c) {
            'M' => .match,
            'I' => .ins,
            'D' => .del,
            'N' => .ref_skip,
            'S' => .soft_clip,
            'H' => .hard_clip,
            'P' => .pad,
            '=' => .equal,
            'X' => .diff,
            'B' => .back,
            else => null,
        };
    }

    pub fn toChar(self: CigarOp) u8 {
        return chars[@intFromEnum(self)];
    }

    /// Whether the operation consumes bases from the query (SEQ): M I S = X.
    pub fn consumesQuery(self: CigarOp) bool {
        return switch (self) {
            .match, .ins, .soft_clip, .equal, .diff => true,
            else => false,
        };
    }

    /// Whether the operation consumes bases from the reference: M D N = X.
    pub fn consumesReference(self: CigarOp) bool {
        return switch (self) {
            .match, .del, .ref_skip, .equal, .diff => true,
            else => false,
        };
    }
};

pub const Cigar = struct { op: CigarOp, len: u32 };

/// Lazily iterates the operations of a CIGAR string, without allocating.
pub const CigarIterator = struct {
    str: []const u8,
    pos: usize = 0,

    pub fn next(self: *CigarIterator) !?Cigar {
        if (self.pos >= self.str.len) return null;
        const start = self.pos;
        var len: u32 = 0;
        while (self.pos < self.str.len and self.str[self.pos] >= '0' and self.str[self.pos] <= '9') : (self.pos += 1) {
            len = std.math.mul(u32, len, 10) catch return error.Overflow;
            len = std.math.add(u32, len, self.str[self.pos] - '0') catch return error.Overflow;
        }
        if (self.pos == start) return error.ParseError; // missing operation length
        if (self.pos >= self.str.len) return error.ParseError; // missing operation char
        const op = CigarOp.fromChar(self.str[self.pos]) orelse return error.ParseError;
        self.pos += 1;
        return .{ .op = op, .len = len };
    }
};

/// Number of query (SEQ) bases spanned by a CIGAR string.
pub fn queryLength(cigar_str: []const u8) !u64 {
    var it = CigarIterator{ .str = cigar_str };
    var sum: u64 = 0;
    while (try it.next()) |c| {
        if (c.op.consumesQuery()) sum += c.len;
    }
    return sum;
}

/// Number of reference bases spanned by a CIGAR string.
pub fn referenceLength(cigar_str: []const u8) !u64 {
    var it = CigarIterator{ .str = cigar_str };
    var sum: u64 = 0;
    while (try it.next()) |c| {
        if (c.op.consumesReference()) sum += c.len;
    }
    return sum;
}

// Optional / auxiliary fields ----------------------------------------------

pub const AuxArray = struct {
    subtype: u8, // one of c C s S i I f
    raw: []const u8, // comma-separated element text

    pub fn elements(self: AuxArray) std.mem.SplitIterator(u8, .scalar) {
        return std.mem.splitScalar(u8, self.raw, ',');
    }
};

pub const AuxValue = union(enum) {
    char: u8, // A
    int: i64, // c C s S i I
    float: f32, // f
    string: []const u8, // Z
    hex: []const u8, // H
    array: AuxArray, // B
};

pub const AuxTag = struct { tag: [2]u8, value: AuxValue };

/// Lazily iterates the TAG:TYPE:VALUE optional fields of an alignment record.
pub const AuxIterator = struct {
    it: std.mem.SplitIterator(u8, .scalar),

    pub fn next(self: *AuxIterator) !?AuxTag {
        while (self.it.next()) |field| {
            if (field.len == 0) continue;
            return try parseAuxField(field);
        }
        return null;
    }
};

fn parseAuxField(field: []const u8) !AuxTag {
    if (field.len < 5 or field[2] != ':' or field[4] != ':') return error.ParseError;
    const tag = [2]u8{ field[0], field[1] };
    const val = field[5..];
    const value: AuxValue = switch (field[3]) {
        'A' => if (val.len == 1) AuxValue{ .char = val[0] } else return error.ParseError,
        'i', 'I', 'c', 'C', 's', 'S' => AuxValue{ .int = std.fmt.parseInt(i64, val, 10) catch return error.ParseError },
        'f' => AuxValue{ .float = std.fmt.parseFloat(f32, val) catch return error.ParseError },
        'Z' => AuxValue{ .string = val },
        'H' => AuxValue{ .hex = val },
        'B' => blk: {
            if (val.len == 0) return error.ParseError;
            const subtype = val[0];
            const raw = if (val.len == 1)
                val[1..] // no elements
            else if (val[1] == ',')
                val[2..]
            else
                return error.ParseError;
            break :blk AuxValue{ .array = .{ .subtype = subtype, .raw = raw } };
        },
        else => return error.ParseError,
    };
    return .{ .tag = tag, .value = value };
}

// Header --------------------------------------------------------------------

pub const Tag = struct { key: [2]u8, value: []const u8 };

/// A @SQ reference sequence entry.
pub const RefSeq = struct {
    name: []const u8,
    len: u64,
};

/// A generic header record (@RG / @PG) retaining all of its tags verbatim.
pub const HeaderLine = struct {
    tags: []Tag,

    pub fn get(self: HeaderLine, key: []const u8) ?[]const u8 {
        for (self.tags) |t| {
            if (std.mem.eql(u8, &t.key, key)) return t.value;
        }
        return null;
    }

    fn free(self: HeaderLine, gpa: Allocator) void {
        for (self.tags) |t| gpa.free(t.value);
        gpa.free(self.tags);
    }
};

pub const Header = struct {
    version_major: u16 = 0,
    version_minor: u16 = 0,
    sort_order: SortOrder = .unknown,
    group_order: GroupAlignment = .none,
    sub_sort: SortOrder = .unknown,
    refs: []RefSeq = &.{},
    read_groups: []HeaderLine = &.{},
    programs: []HeaderLine = &.{},
    comments: [][]const u8 = &.{},
    /// Reference name -> tid (0-based index into `refs`).
    name2tid: std.StringHashMapUnmanaged(i32) = .empty,

    /// Resolve a reference name to its tid, or null if not in the dictionary.
    pub fn tid(self: *const Header, name: []const u8) ?i32 {
        return self.name2tid.get(name);
    }

    pub fn refName(self: *const Header, t: i32) ?[]const u8 {
        if (t < 0) return null;
        const idx: usize = @intCast(t);
        if (idx >= self.refs.len) return null;
        return self.refs[idx].name;
    }

    pub fn refLen(self: *const Header, t: i32) ?u64 {
        if (t < 0) return null;
        const idx: usize = @intCast(t);
        if (idx >= self.refs.len) return null;
        return self.refs[idx].len;
    }

    pub fn deinit(self: *Header, gpa: Allocator) void {
        self.name2tid.deinit(gpa);
        for (self.refs) |rs| gpa.free(rs.name);
        gpa.free(self.refs);
        for (self.read_groups) |hl| hl.free(gpa);
        gpa.free(self.read_groups);
        for (self.programs) |hl| hl.free(gpa);
        gpa.free(self.programs);
        for (self.comments) |c| gpa.free(c);
        gpa.free(self.comments);
        self.* = .{};
    }
};

// Alignment record ----------------------------------------------------------

/// A parsed alignment record. Slice fields point into the reader's line buffer
/// and are only valid until the next `SamReader.next` call.
pub const Alignment = struct {
    qname: []const u8,
    flag: u16,
    rname: []const u8, // raw RNAME text ("*" preserved)
    tid: i32, // resolved reference id, -1 if unmapped/unknown
    pos: i64, // 1-based leftmost position as in the file; 0 = none
    mapq: u8,
    cigar_str: []const u8, // empty when CIGAR is "*"
    rnext: []const u8, // raw RNEXT text ("=" / "*" preserved)
    mtid: i32, // resolved mate reference id, -1 if unmapped/unknown
    pnext: i64, // 1-based mate position; 0 = none
    tlen: i64, // observed template length
    seq: []const u8, // empty when SEQ is "*"
    qual: []const u8, // raw QUAL (Phred+33); empty when "*"
    aux_str: []const u8, // remaining tab-separated optional fields

    pub fn flags(self: Alignment) SamAlignmentFlags {
        return SamAlignmentFlags.toSamAlignmentFlags(self.flag);
    }

    pub fn cigar(self: Alignment) CigarIterator {
        return .{ .str = self.cigar_str };
    }

    pub fn aux(self: Alignment) AuxIterator {
        return .{ .it = std.mem.splitScalar(u8, self.aux_str, '\t') };
    }

    /// Phred quality of base `i` (QUAL byte minus 33). Asserts QUAL is present.
    pub fn qualPhred(self: Alignment, i: usize) u8 {
        return self.qual[i] - 33;
    }
};

// Reader --------------------------------------------------------------------

pub const SamReader = struct {
    allocator: Allocator,
    r: *std.Io.Reader,

    /// Parse all leading `@` header lines into an owned `Header`, stopping at
    /// the first alignment line (or end of stream). Caller owns the returned
    /// header and must `deinit` it.
    pub fn readHeader(self: *SamReader) !Header {
        const gpa = self.allocator;

        var meta: Header = .{};

        var refs: std.ArrayListUnmanaged(RefSeq) = .empty;
        errdefer {
            for (refs.items) |rs| gpa.free(rs.name);
            refs.deinit(gpa);
        }
        var name2tid: std.StringHashMapUnmanaged(i32) = .empty;
        errdefer name2tid.deinit(gpa);
        var read_groups: std.ArrayListUnmanaged(HeaderLine) = .empty;
        errdefer {
            for (read_groups.items) |hl| hl.free(gpa);
            read_groups.deinit(gpa);
        }
        var programs: std.ArrayListUnmanaged(HeaderLine) = .empty;
        errdefer {
            for (programs.items) |hl| hl.free(gpa);
            programs.deinit(gpa);
        }
        var comments: std.ArrayListUnmanaged([]const u8) = .empty;
        errdefer {
            for (comments.items) |c| gpa.free(c);
            comments.deinit(gpa);
        }

        while (true) {
            const b = self.r.peekByte() catch |err| switch (err) {
                error.EndOfStream => break,
                else => |e| return e,
            };
            if (b != '@') break; // start of the alignment section

            const raw = (try self.r.takeDelimiter('\n')) orelse break;
            const line = std.mem.trimEnd(u8, raw, "\r");
            if (line.len < 3) continue;
            const rec = line[0..3];

            if (std.mem.eql(u8, rec, "@HD")) {
                try parseHd(line, &meta);
            } else if (std.mem.eql(u8, rec, "@SQ")) {
                var name: ?[]const u8 = null;
                var ln: u64 = 0;
                var cols = std.mem.splitScalar(u8, line, '\t');
                _ = cols.next(); // @SQ
                while (cols.next()) |col| {
                    if (col.len < 3 or col[2] != ':') continue;
                    const key = col[0..2];
                    const val = col[3..];
                    if (std.mem.eql(u8, key, "SN")) {
                        name = val;
                    } else if (std.mem.eql(u8, key, "LN")) {
                        ln = std.fmt.parseUnsigned(u64, val, 10) catch return error.ParseError;
                    }
                }
                const sn = name orelse return error.ParseError;
                const owned = try gpa.dupe(u8, sn);
                {
                    errdefer gpa.free(owned);
                    try refs.append(gpa, .{ .name = owned, .len = ln });
                }
                // `owned` is now owned by `refs`; its errdefer covers it below.
                const t: i32 = @intCast(refs.items.len - 1);
                try name2tid.put(gpa, owned, t);
            } else if (std.mem.eql(u8, rec, "@RG")) {
                const hl = try parseHeaderLine(gpa, line);
                errdefer hl.free(gpa);
                try read_groups.append(gpa, hl);
            } else if (std.mem.eql(u8, rec, "@PG")) {
                const hl = try parseHeaderLine(gpa, line);
                errdefer hl.free(gpa);
                try programs.append(gpa, hl);
            } else if (std.mem.eql(u8, rec, "@CO")) {
                const comment = if (line.len > 4) line[4..] else "";
                const owned = try gpa.dupe(u8, comment);
                {
                    errdefer gpa.free(owned);
                    try comments.append(gpa, owned);
                }
            }
            // Unknown @ records are ignored.
        }

        meta.refs = try refs.toOwnedSlice(gpa);
        meta.read_groups = try read_groups.toOwnedSlice(gpa);
        meta.programs = try programs.toOwnedSlice(gpa);
        meta.comments = try comments.toOwnedSlice(gpa);
        meta.name2tid = name2tid;
        return meta;
    }

    /// Parse the next alignment record, or null at end of stream. The returned
    /// `Alignment` borrows the reader's line buffer until the next call.
    pub fn next(self: *SamReader, header: *const Header) !?Alignment {
        while (true) {
            const raw = (try self.r.takeDelimiter('\n')) orelse return null;
            const line = std.mem.trimEnd(u8, raw, "\r");
            if (line.len == 0) continue; // tolerate blank lines
            return try parseAlignment(header, line);
        }
    }
};

pub fn reader(allocator: Allocator, r: *std.Io.Reader) SamReader {
    return .{ .allocator = allocator, .r = r };
}

fn parseSortOrder(s: []const u8) ?SortOrder {
    if (std.mem.eql(u8, s, "unknown")) return .unknown;
    if (std.mem.eql(u8, s, "unsorted")) return .unsorted;
    if (std.mem.eql(u8, s, "queryname")) return .query_name;
    if (std.mem.eql(u8, s, "coordinate")) return .coordinate;
    return null;
}

fn parseGroupOrder(s: []const u8) ?GroupAlignment {
    if (std.mem.eql(u8, s, "none")) return .none;
    if (std.mem.eql(u8, s, "query")) return .query;
    if (std.mem.eql(u8, s, "reference")) return .reference;
    return null;
}

fn parseHd(line: []const u8, meta: *Header) !void {
    var cols = std.mem.splitScalar(u8, line, '\t');
    _ = cols.next(); // @HD
    while (cols.next()) |col| {
        if (col.len < 3 or col[2] != ':') continue;
        const key = col[0..2];
        const val = col[3..]; // value may itself contain ':'
        if (std.mem.eql(u8, key, "VN")) {
            var vi = std.mem.splitScalar(u8, val, '.');
            const maj = vi.next() orelse return error.ParseError;
            meta.version_major = std.fmt.parseUnsigned(u16, maj, 10) catch return error.ParseError;
            if (vi.next()) |min| {
                meta.version_minor = std.fmt.parseUnsigned(u16, min, 10) catch return error.ParseError;
            }
        } else if (std.mem.eql(u8, key, "SO")) {
            meta.sort_order = parseSortOrder(val) orelse return error.ParseError;
        } else if (std.mem.eql(u8, key, "GO")) {
            meta.group_order = parseGroupOrder(val) orelse return error.ParseError;
        } else if (std.mem.eql(u8, key, "SS")) {
            // SS is "sort-order(:sub-sort)+"; classify the leading component.
            var si = std.mem.splitScalar(u8, val, ':');
            const first = si.next() orelse return error.ParseError;
            meta.sub_sort = parseSortOrder(first) orelse return error.ParseError;
        }
    }
}

fn parseHeaderLine(gpa: Allocator, line: []const u8) !HeaderLine {
    var tags: std.ArrayListUnmanaged(Tag) = .empty;
    errdefer {
        for (tags.items) |t| gpa.free(t.value);
        tags.deinit(gpa);
    }
    var cols = std.mem.splitScalar(u8, line, '\t');
    _ = cols.next(); // record type (@RG / @PG)
    while (cols.next()) |col| {
        if (col.len < 3 or col[2] != ':') continue;
        const key = [2]u8{ col[0], col[1] };
        const val = try gpa.dupe(u8, col[3..]);
        {
            errdefer gpa.free(val);
            try tags.append(gpa, .{ .key = key, .value = val });
        }
    }
    return .{ .tags = try tags.toOwnedSlice(gpa) };
}

fn parseAlignment(header: *const Header, line: []const u8) !Alignment {
    var it = std.mem.splitScalar(u8, line, '\t');
    const qname = it.next() orelse return error.ParseError;
    const flag_s = it.next() orelse return error.ParseError;
    const rname = it.next() orelse return error.ParseError;
    const pos_s = it.next() orelse return error.ParseError;
    const mapq_s = it.next() orelse return error.ParseError;
    const cigar_s = it.next() orelse return error.ParseError;
    const rnext = it.next() orelse return error.ParseError;
    const pnext_s = it.next() orelse return error.ParseError;
    const tlen_s = it.next() orelse return error.ParseError;
    const seq_s = it.next() orelse return error.ParseError;
    const qual_s = it.next() orelse return error.ParseError;
    const aux_str = it.rest(); // remaining optional fields (may be empty)

    const flag = std.fmt.parseUnsigned(u16, flag_s, 10) catch return error.ParseError;
    const mapq = std.fmt.parseUnsigned(u8, mapq_s, 10) catch return error.ParseError;
    const pos = std.fmt.parseInt(i64, pos_s, 10) catch return error.ParseError;
    const pnext = std.fmt.parseInt(i64, pnext_s, 10) catch return error.ParseError;
    const tlen = std.fmt.parseInt(i64, tlen_s, 10) catch return error.ParseError;

    const tid: i32 = if (std.mem.eql(u8, rname, "*")) -1 else (header.tid(rname) orelse -1);
    const mtid: i32 = if (std.mem.eql(u8, rnext, "="))
        tid
    else if (std.mem.eql(u8, rnext, "*"))
        -1
    else
        (header.tid(rnext) orelse -1);

    const cigar_str: []const u8 = if (std.mem.eql(u8, cigar_s, "*")) "" else cigar_s;
    const seq: []const u8 = if (std.mem.eql(u8, seq_s, "*")) "" else seq_s;
    const qual: []const u8 = if (std.mem.eql(u8, qual_s, "*")) "" else qual_s;

    // Consistency checks: CIGAR query span must equal SEQ length (when both are
    // present), and QUAL length must equal SEQ length (when both are present).
    if (seq.len != 0 and cigar_str.len != 0) {
        if ((try queryLength(cigar_str)) != seq.len) return error.ParseError;
    }
    if (seq.len != 0 and qual.len != 0 and qual.len != seq.len) return error.ParseError;

    return .{
        .qname = qname,
        .flag = flag,
        .rname = rname,
        .tid = tid,
        .pos = pos,
        .mapq = mapq,
        .cigar_str = cigar_str,
        .rnext = rnext,
        .mtid = mtid,
        .pnext = pnext,
        .tlen = tlen,
        .seq = seq,
        .qual = qual,
        .aux_str = aux_str,
    };
}

// Tests ---------------------------------------------------------------------

test "flags round-trip" {
    const f = SamAlignmentFlags.toSamAlignmentFlags(99); // paired|proper|mate_reverse|read1
    try testing.expect(f.paired and f.proper_pair and f.mate_reverse and f.read1);
    try testing.expect(!f.reverse and !f.unmap);
    try testing.expectEqual(@as(u16, 99), f.toValue());
}

test "parse header: HD/SQ/RG/PG/CO and reference dictionary" {
    const text =
        "@HD\tVN:1.6\tSO:coordinate\n" ++
        "@SQ\tSN:chr1\tLN:1000\n" ++
        "@SQ\tSN:chr2\tLN:2000\n" ++
        "@RG\tID:rg1\tSM:sample1\tPL:ILLUMINA\n" ++
        "@PG\tID:bwa\tPN:bwa\tVN:0.7.17\n" ++
        "@CO\tthis is a comment\n";
    var r: std.Io.Reader = .fixed(text);
    var sr = reader(testing.allocator, &r);
    var header = try sr.readHeader();
    defer header.deinit(testing.allocator);

    try testing.expectEqual(@as(u16, 1), header.version_major);
    try testing.expectEqual(@as(u16, 6), header.version_minor);
    try testing.expectEqual(SortOrder.coordinate, header.sort_order);

    try testing.expectEqual(@as(usize, 2), header.refs.len);
    try testing.expectEqual(@as(?i32, 0), header.tid("chr1"));
    try testing.expectEqual(@as(?i32, 1), header.tid("chr2"));
    try testing.expectEqual(@as(?i32, null), header.tid("chrX"));
    try testing.expectEqualStrings("chr1", header.refName(0).?);
    try testing.expectEqual(@as(?u64, 2000), header.refLen(1));

    try testing.expectEqual(@as(usize, 1), header.read_groups.len);
    try testing.expectEqualStrings("sample1", header.read_groups[0].get("SM").?);
    try testing.expectEqual(@as(?[]const u8, null), header.read_groups[0].get("XX"));
    try testing.expectEqualStrings("bwa", header.programs[0].get("ID").?);

    try testing.expectEqual(@as(usize, 1), header.comments.len);
    try testing.expectEqualStrings("this is a comment", header.comments[0]);
}

test "parse alignment: fields, cigar, flags, aux, sentinels" {
    const text =
        "@HD\tVN:1.6\tSO:coordinate\n" ++
        "@SQ\tSN:chr1\tLN:1000\n" ++
        "@SQ\tSN:chr2\tLN:2000\n" ++
        "r1\t99\tchr1\t100\t60\t8M2I5M\t=\t200\t150\tACGTACGTACGTACG\tIIIIIIIIIIIIIII\tNM:i:2\tMD:Z:15\n" ++
        "r2\t4\t*\t0\t0\t*\t*\t0\t0\t*\t*\n";
    var r: std.Io.Reader = .fixed(text);
    var sr = reader(testing.allocator, &r);
    var header = try sr.readHeader();
    defer header.deinit(testing.allocator);

    const a1 = (try sr.next(&header)).?;
    try testing.expectEqualStrings("r1", a1.qname);
    try testing.expectEqual(@as(u16, 99), a1.flag);
    try testing.expectEqual(@as(i32, 0), a1.tid); // chr1
    try testing.expectEqual(@as(i64, 100), a1.pos);
    try testing.expectEqual(@as(u8, 60), a1.mapq);
    try testing.expectEqual(@as(i32, 0), a1.mtid); // "=" -> same as chr1
    try testing.expectEqual(@as(i64, 200), a1.pnext);
    try testing.expectEqual(@as(i64, 150), a1.tlen);
    try testing.expectEqualStrings("ACGTACGTACGTACG", a1.seq);

    const f = a1.flags();
    try testing.expect(f.paired and f.proper_pair and f.mate_reverse and f.read1);
    try testing.expect(!f.unmap);

    var ci = a1.cigar();
    const c0 = (try ci.next()).?;
    try testing.expectEqual(CigarOp.match, c0.op);
    try testing.expectEqual(@as(u32, 8), c0.len);
    const c1 = (try ci.next()).?;
    try testing.expectEqual(CigarOp.ins, c1.op);
    try testing.expectEqual(@as(u32, 2), c1.len);
    try testing.expectEqual(@as(u64, 15), try queryLength(a1.cigar_str));
    try testing.expectEqual(@as(u64, 13), try referenceLength(a1.cigar_str));

    var ai = a1.aux();
    const nm = (try ai.next()).?;
    try testing.expectEqualStrings("NM", &nm.tag);
    try testing.expectEqual(@as(i64, 2), nm.value.int);
    const md = (try ai.next()).?;
    try testing.expectEqualStrings("MD", &md.tag);
    try testing.expectEqualStrings("15", md.value.string);
    try testing.expect((try ai.next()) == null);

    const a2 = (try sr.next(&header)).?;
    try testing.expectEqual(@as(i32, -1), a2.tid);
    try testing.expect(a2.flags().unmap);
    try testing.expectEqual(@as(usize, 0), a2.cigar_str.len);
    try testing.expectEqual(@as(usize, 0), a2.seq.len);
    try testing.expectEqual(@as(usize, 0), a2.qual.len);

    try testing.expect((try sr.next(&header)) == null);
}

test "aux types A, f, and B array" {
    const text =
        "@SQ\tSN:c\tLN:100\n" ++
        "x\t0\tc\t1\t0\t4M\t*\t0\t0\tACGT\tIIII\tXA:A:P\tXF:f:3.5\tXB:B:c,1,-2,3\n";
    var r: std.Io.Reader = .fixed(text);
    var sr = reader(testing.allocator, &r);
    var header = try sr.readHeader();
    defer header.deinit(testing.allocator);

    const a = (try sr.next(&header)).?;
    var ai = a.aux();
    const xa = (try ai.next()).?;
    try testing.expectEqual(@as(u8, 'P'), xa.value.char);
    const xf = (try ai.next()).?;
    try testing.expectApproxEqAbs(@as(f32, 3.5), xf.value.float, 0.0001);
    const xb = (try ai.next()).?;
    try testing.expectEqual(@as(u8, 'c'), xb.value.array.subtype);
    var elems = xb.value.array.elements();
    try testing.expectEqualStrings("1", elems.next().?);
    try testing.expectEqualStrings("-2", elems.next().?);
    try testing.expectEqualStrings("3", elems.next().?);
    try testing.expect(elems.next() == null);
}

test "CIGAR/SEQ length mismatch is a parse error" {
    const text =
        "@SQ\tSN:c\tLN:100\n" ++
        "x\t0\tc\t1\t0\t10M\t*\t0\t0\tACGT\tIIII\n"; // CIGAR spans 10 query bases, SEQ has 4
    var r: std.Io.Reader = .fixed(text);
    var sr = reader(testing.allocator, &r);
    var header = try sr.readHeader();
    defer header.deinit(testing.allocator);
    try testing.expectError(error.ParseError, sr.next(&header));
}

test "read embedded t1.sam end to end" {
    const file = @embedFile("./test/t1.sam");
    var r: std.Io.Reader = .fixed(file);
    var sr = reader(testing.allocator, &r);
    var header = try sr.readHeader();
    defer header.deinit(testing.allocator);

    try testing.expectEqual(@as(usize, 2), header.refs.len);
    try testing.expectEqual(@as(?i32, 0), header.tid("ref1"));
    try testing.expectEqual(@as(?i32, 1), header.tid("ref2"));

    var count: usize = 0;
    while (try sr.next(&header)) |rec| {
        // Every mapped record's CIGAR query span already validated in next();
        // here just exercise the reference dictionary resolution.
        if (rec.tid >= 0) {
            try testing.expect(header.refName(rec.tid) != null);
        }
        count += 1;
    }
    try testing.expectEqual(@as(usize, 5), count);
}
