const std = @import("std");
const sam = @import("sam.zig");
const testing = std.testing;
const Allocator = std.mem.Allocator;
const flate = std.compress.flate;

// BAM: the binary, BGZF-compressed encoding of SAM. This module provides a
// reader (decompress + parse) and a writer (encode + compress), reusing the
// shared value types in sam.zig (SamAlignmentFlags, CigarOp, Header, ...).
//
// Reference: https://samtools.github.io/hts-specs/SAMv1.pdf (§4 BAM, §4.1 BGZF)

pub const magic = "BAM\x01";

/// 4-bit SEQ code -> base character (BAM nibble encoding).
const seq_nt16_str = "=ACMGRSVTWYHKDBN";

/// Base character -> 4-bit SEQ code (inverse of seq_nt16_str, case-insensitive).
fn seqNt16(base: u8) u4 {
    return switch (base) {
        '=' => 0,
        'A', 'a' => 1,
        'C', 'c' => 2,
        'M', 'm' => 3,
        'G', 'g' => 4,
        'R', 'r' => 5,
        'S', 's' => 6,
        'V', 'v' => 7,
        'T', 't' => 8,
        'W', 'w' => 9,
        'Y', 'y' => 10,
        'H', 'h' => 11,
        'K', 'k' => 12,
        'D', 'd' => 13,
        'B', 'b' => 14,
        else => 15, // N / unknown
    };
}

// BGZF ----------------------------------------------------------------------

/// Standard 28-byte BGZF end-of-file marker (an empty block).
pub const bgzf_eof = [_]u8{
    0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00,
    0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00,
    0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00,
};

/// Maximum uncompressed payload per BGZF block (htslib's 0xff00).
const max_block_payload = 0xff00;

/// Decompress an entire BGZF stream into a single owned buffer. BGZF is a
/// series of concatenated gzip members; `flate.Decompress` only decodes one
/// member per call, so each block is framed and inflated individually here.
pub fn inflateAll(gpa: Allocator, input: *std.Io.Reader) ![]u8 {
    var out: std.ArrayList(u8) = .empty;
    errdefer out.deinit(gpa);

    while (true) {
        const hdr = input.takeArray(12) catch |err| switch (err) {
            error.EndOfStream => break, // clean end of the BGZF stream
            else => |e| return e,
        };
        if (hdr[0] != 0x1f or hdr[1] != 0x8b or hdr[2] != 0x08 or (hdr[3] & 0x04) == 0)
            return error.InvalidBgzf;
        const xlen = std.mem.readInt(u16, hdr[10..12], .little);

        const extra = try input.take(xlen);
        var bsize: ?u16 = null;
        var i: usize = 0;
        while (i + 4 <= extra.len) {
            const slen = std.mem.readInt(u16, extra[i + 2 ..][0..2], .little);
            if (extra[i] == 'B' and extra[i + 1] == 'C' and slen == 2) {
                if (i + 6 > extra.len) return error.InvalidBgzf;
                bsize = std.mem.readInt(u16, extra[i + 4 ..][0..2], .little);
            }
            i += 4 + slen;
        }
        const bs = bsize orelse return error.InvalidBgzf;

        const total: usize = @as(usize, bs) + 1;
        if (total < 12 + @as(usize, xlen) + 8) return error.InvalidBgzf;
        const cdata_len = total - 12 - @as(usize, xlen) - 8;

        const cdata = try input.take(cdata_len);
        var cin: std.Io.Reader = .fixed(cdata);
        var dec = flate.Decompress.init(&cin, .raw, &.{});
        try dec.reader.appendRemaining(gpa, &out, .unlimited);

        _ = try input.takeArray(8); // skip CRC32 + ISIZE footer
    }

    return out.toOwnedSlice(gpa);
}

/// Compress `payload` (<= max_block_payload bytes) as one BGZF block to `out`.
pub fn writeBlock(out: *std.Io.Writer, payload: []const u8) !void {
    var cbuf: [0x10000]u8 = undefined;
    var cw: std.Io.Writer = .fixed(&cbuf);
    var window: [flate.max_window_len]u8 = undefined;
    var comp = try flate.Compress.init(&cw, &window, .raw, .default);
    try comp.writer.writeAll(payload);
    try comp.finish();
    const cdata = cw.buffered();

    const total = 18 + cdata.len + 8;
    var hdr = [18]u8{
        0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00,
        0x00, 0xff, 0x06, 0x00, 'B',  'C',  0x02, 0x00,
        0x00, 0x00,
    };
    std.mem.writeInt(u16, hdr[16..18], @intCast(total - 1), .little);
    try out.writeAll(&hdr);
    try out.writeAll(cdata);

    var footer: [8]u8 = undefined;
    std.mem.writeInt(u32, footer[0..4], std.hash.Crc32.hash(payload), .little);
    std.mem.writeInt(u32, footer[4..8], @intCast(payload.len), .little);
    try out.writeAll(&footer);
}

/// Emit the standard BGZF end-of-file marker block.
pub fn writeEof(out: *std.Io.Writer) !void {
    try out.writeAll(&bgzf_eof);
}

// Cursor over the decompressed BAM byte stream -------------------------------

const Cursor = struct {
    buf: []const u8,
    pos: usize = 0,

    fn takeBytes(self: *Cursor, n: usize) ![]const u8 {
        if (self.pos + n > self.buf.len) return error.ParseError;
        const s = self.buf[self.pos..][0..n];
        self.pos += n;
        return s;
    }
    fn byte(self: *Cursor) !u8 {
        return (try self.takeBytes(1))[0];
    }
    fn u16le(self: *Cursor) !u16 {
        return std.mem.readInt(u16, (try self.takeBytes(2))[0..2], .little);
    }
    fn u32le(self: *Cursor) !u32 {
        return std.mem.readInt(u32, (try self.takeBytes(4))[0..4], .little);
    }
    fn i32le(self: *Cursor) !i32 {
        return std.mem.readInt(i32, (try self.takeBytes(4))[0..4], .little);
    }
};

// Alignment record ----------------------------------------------------------

/// CIGAR iterator over the raw BAM cigar bytes (n_cigar_op little-endian u32).
pub const CigarIterator = struct {
    raw: []const u8,
    pos: usize = 0,

    pub fn next(self: *CigarIterator) ?sam.Cigar {
        if (self.pos + 4 > self.raw.len) return null;
        const v = std.mem.readInt(u32, self.raw[self.pos..][0..4], .little);
        self.pos += 4;
        return .{ .op = @enumFromInt(@as(u8, @intCast(v & 0xf))), .len = v >> 4 };
    }
};

/// A binary B-array optional value (subtype + count + raw little-endian elements).
pub const AuxArray = struct {
    subtype: u8, // c C s S i I f
    count: u32,
    raw: []const u8,

    pub fn elemSize(self: AuxArray) usize {
        return switch (self.subtype) {
            'c', 'C' => 1,
            's', 'S' => 2,
            'i', 'I', 'f' => 4,
            else => 0,
        };
    }
    pub fn intAt(self: AuxArray, i: usize) i64 {
        const off = i * self.elemSize();
        return switch (self.subtype) {
            'c' => @as(i8, @bitCast(self.raw[off])),
            'C' => self.raw[off],
            's' => std.mem.readInt(i16, self.raw[off..][0..2], .little),
            'S' => std.mem.readInt(u16, self.raw[off..][0..2], .little),
            'i' => std.mem.readInt(i32, self.raw[off..][0..4], .little),
            'I' => std.mem.readInt(u32, self.raw[off..][0..4], .little),
            else => 0,
        };
    }
    pub fn floatAt(self: AuxArray, i: usize) f32 {
        return @bitCast(std.mem.readInt(u32, self.raw[i * 4 ..][0..4], .little));
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

/// Iterates the binary optional fields (tag + type + value) of a record.
pub const AuxIterator = struct {
    raw: []const u8,
    pos: usize = 0,

    pub fn next(self: *AuxIterator) !?AuxTag {
        if (self.pos >= self.raw.len) return null;
        if (self.pos + 3 > self.raw.len) return error.ParseError;
        const tag = [2]u8{ self.raw[self.pos], self.raw[self.pos + 1] };
        const typ = self.raw[self.pos + 2];
        self.pos += 3;
        const value: AuxValue = switch (typ) {
            'A' => .{ .char = try self.readByte() },
            'c' => .{ .int = @as(i8, @bitCast(try self.readByte())) },
            'C' => .{ .int = try self.readByte() },
            's' => .{ .int = try self.readIntT(i16) },
            'S' => .{ .int = try self.readIntT(u16) },
            'i' => .{ .int = try self.readIntT(i32) },
            'I' => .{ .int = try self.readIntT(u32) },
            'f' => .{ .float = @bitCast(try self.readIntT(u32)) },
            'Z' => .{ .string = try self.readCStr() },
            'H' => .{ .hex = try self.readCStr() },
            'B' => try self.readArray(),
            else => return error.ParseError,
        };
        return .{ .tag = tag, .value = value };
    }

    fn readByte(self: *AuxIterator) !u8 {
        if (self.pos >= self.raw.len) return error.ParseError;
        defer self.pos += 1;
        return self.raw[self.pos];
    }
    fn readIntT(self: *AuxIterator, comptime T: type) !T {
        const n = @sizeOf(T);
        if (self.pos + n > self.raw.len) return error.ParseError;
        const v = std.mem.readInt(T, self.raw[self.pos..][0..n], .little);
        self.pos += n;
        return v;
    }
    fn readCStr(self: *AuxIterator) ![]const u8 {
        const start = self.pos;
        while (self.pos < self.raw.len and self.raw[self.pos] != 0) self.pos += 1;
        if (self.pos >= self.raw.len) return error.ParseError; // missing NUL
        const s = self.raw[start..self.pos];
        self.pos += 1;
        return s;
    }
    fn readArray(self: *AuxIterator) !AuxValue {
        const subtype = try self.readByte();
        const count = try self.readIntT(u32);
        const es: usize = switch (subtype) {
            'c', 'C' => 1,
            's', 'S' => 2,
            'i', 'I', 'f' => 4,
            else => return error.ParseError,
        };
        const nbytes = @as(usize, count) * es;
        if (self.pos + nbytes > self.raw.len) return error.ParseError;
        const raw = self.raw[self.pos..][0..nbytes];
        self.pos += nbytes;
        return .{ .array = .{ .subtype = subtype, .count = count, .raw = raw } };
    }
};

/// A parsed BAM alignment record. Slice fields point into the reader's owned
/// decompressed buffer and stay valid for the reader's lifetime.
pub const Record = struct {
    ref_id: i32, // -1 if unmapped/unknown
    pos: i64, // 1-based (BAM 0-based +1; unmapped -1 -> 0), matching the SAM reader
    mapq: u8,
    bin: u16,
    flag: u16,
    next_ref_id: i32,
    next_pos: i64, // 1-based
    tlen: i64,
    l_seq: u32,
    qname: []const u8, // without the trailing NUL
    cigar_raw: []const u8, // n_cigar_op * 4 bytes
    seq_packed: []const u8, // (l_seq + 1) / 2 bytes, 4-bit packed
    qual: []const u8, // l_seq bytes, raw Phred (0xff = unavailable)
    aux_raw: []const u8,

    pub fn flags(self: Record) sam.SamAlignmentFlags {
        return sam.SamAlignmentFlags.toSamAlignmentFlags(self.flag);
    }
    pub fn cigar(self: Record) CigarIterator {
        return .{ .raw = self.cigar_raw };
    }
    pub fn aux(self: Record) AuxIterator {
        return .{ .raw = self.aux_raw };
    }
    /// Decode base `i` of SEQ to its character.
    pub fn seqBase(self: Record, i: usize) u8 {
        const b = self.seq_packed[i >> 1];
        const nib: u4 = if (i & 1 == 0) @intCast(b >> 4) else @intCast(b & 0xf);
        return seq_nt16_str[nib];
    }
    /// True when QUAL is unavailable (BAM stores 0xff in every byte).
    pub fn qualUnavailable(self: Record) bool {
        return self.qual.len > 0 and self.qual[0] == 0xff;
    }
};

fn parseRecord(cur: *Cursor) !?Record {
    if (cur.pos >= cur.buf.len) return null;
    const block_size = try cur.u32le();
    const start = cur.pos;

    const ref_id = try cur.i32le();
    const pos0 = try cur.i32le();
    const l_read_name = try cur.byte();
    const mapq = try cur.byte();
    const bin = try cur.u16le();
    const n_cigar_op = try cur.u16le();
    const flag = try cur.u16le();
    const l_seq = try cur.u32le();
    const next_ref_id = try cur.i32le();
    const next_pos0 = try cur.i32le();
    const tlen = try cur.i32le();

    if (l_read_name == 0) return error.ParseError;
    const qname_z = try cur.takeBytes(l_read_name);
    const qname = qname_z[0 .. l_read_name - 1];
    const cigar_raw = try cur.takeBytes(@as(usize, n_cigar_op) * 4);
    const seq_packed = try cur.takeBytes((@as(usize, l_seq) + 1) / 2);
    const qual = try cur.takeBytes(l_seq);

    const consumed = cur.pos - start;
    if (consumed > block_size) return error.ParseError;
    const aux_raw = try cur.takeBytes(block_size - consumed);

    return .{
        .ref_id = ref_id,
        .pos = if (pos0 < 0) 0 else @as(i64, pos0) + 1,
        .mapq = mapq,
        .bin = bin,
        .flag = flag,
        .next_ref_id = next_ref_id,
        .next_pos = if (next_pos0 < 0) 0 else @as(i64, next_pos0) + 1,
        .tlen = tlen,
        .l_seq = l_seq,
        .qname = qname,
        .cigar_raw = cigar_raw,
        .seq_packed = seq_packed,
        .qual = qual,
        .aux_raw = aux_raw,
    };
}

// Reader --------------------------------------------------------------------

pub const BamReader = struct {
    gpa: Allocator,
    decompressed: []u8, // owned
    header: sam.Header, // owned
    header_text: []const u8, // SAM header text, slice into `decompressed`
    pos: usize, // cursor at the next record

    pub fn next(self: *BamReader) !?Record {
        var cur = Cursor{ .buf = self.decompressed, .pos = self.pos };
        const rec = (try parseRecord(&cur)) orelse return null;
        self.pos = cur.pos;
        return rec;
    }

    pub fn deinit(self: *BamReader) void {
        self.header.deinit(self.gpa);
        self.gpa.free(self.decompressed);
        self.* = undefined;
    }
};

/// Open a BAM stream: decompress it, parse the magic/header/reference list, and
/// position at the first alignment record. Caller owns the returned reader.
pub fn open(gpa: Allocator, input: *std.Io.Reader) !BamReader {
    const decompressed = try inflateAll(gpa, input);
    errdefer gpa.free(decompressed);

    var cur = Cursor{ .buf = decompressed };
    if (!std.mem.eql(u8, try cur.takeBytes(4), magic)) return error.NotBam;
    const l_text = try cur.u32le();
    const text = try cur.takeBytes(l_text);
    const n_ref = try cur.u32le();

    // Parse @HD/@RG/@PG/@CO (and @SQ) from the header text via the SAM parser...
    var text_reader: std.Io.Reader = .fixed(text);
    var shr = sam.reader(gpa, &text_reader);
    var header = try shr.readHeader();
    errdefer header.deinit(gpa);

    // ...but the binary reference list is authoritative for tid order, so
    // rebuild the reference dictionary from it.
    for (header.refs) |rs| gpa.free(rs.name);
    gpa.free(header.refs);
    header.refs = &.{};
    header.name2tid.deinit(gpa);
    header.name2tid = .empty;

    var refs: std.ArrayList(sam.RefSeq) = .empty;
    errdefer {
        for (refs.items) |rs| gpa.free(rs.name);
        refs.deinit(gpa);
    }
    var name2tid: std.StringHashMapUnmanaged(i32) = .empty;
    errdefer name2tid.deinit(gpa);

    var ri: usize = 0;
    while (ri < n_ref) : (ri += 1) {
        const l_name = try cur.u32le();
        const name_z = try cur.takeBytes(l_name);
        const name = if (l_name > 0) name_z[0 .. l_name - 1] else name_z;
        const l_ref = try cur.u32le();
        const owned = try gpa.dupe(u8, name);
        {
            errdefer gpa.free(owned);
            try refs.append(gpa, .{ .name = owned, .len = l_ref });
        }
        try name2tid.put(gpa, owned, @intCast(ri));
    }
    header.refs = try refs.toOwnedSlice(gpa);
    header.name2tid = name2tid;

    return .{
        .gpa = gpa,
        .decompressed = decompressed,
        .header = header,
        .header_text = text,
        .pos = cur.pos,
    };
}

// Writer --------------------------------------------------------------------

pub const Writer = struct {
    gpa: Allocator,
    out: *std.Io.Writer,
    buf: std.ArrayList(u8), // staged uncompressed BAM bytes

    pub fn writeHeader(self: *Writer, text: []const u8, refs: []const sam.RefSeq) !void {
        try self.append(magic);
        try self.appendU32(@intCast(text.len));
        try self.append(text);
        try self.appendU32(@intCast(refs.len));
        for (refs) |rs| {
            try self.appendU32(@intCast(rs.name.len + 1));
            try self.append(rs.name);
            try self.appendByte(0);
            try self.appendU32(@intCast(rs.len));
        }
        try self.flushBlocks(false);
    }

    pub fn writeRecord(self: *Writer, rec: *const Record) !void {
        const l_read_name: u8 = @intCast(rec.qname.len + 1);
        const n_cigar_op: u16 = @intCast(rec.cigar_raw.len / 4);
        const block_size: u32 = @intCast(32 + @as(usize, l_read_name) +
            rec.cigar_raw.len + rec.seq_packed.len + rec.qual.len + rec.aux_raw.len);

        try self.appendU32(block_size);
        try self.appendI32(rec.ref_id);
        try self.appendI32(if (rec.pos == 0) -1 else @intCast(rec.pos - 1));
        try self.appendByte(l_read_name);
        try self.appendByte(rec.mapq);
        try self.appendU16(rec.bin);
        try self.appendU16(n_cigar_op);
        try self.appendU16(rec.flag);
        try self.appendU32(rec.l_seq);
        try self.appendI32(rec.next_ref_id);
        try self.appendI32(if (rec.next_pos == 0) -1 else @intCast(rec.next_pos - 1));
        try self.appendI32(@intCast(rec.tlen));
        try self.append(rec.qname);
        try self.appendByte(0);
        try self.append(rec.cigar_raw);
        try self.append(rec.seq_packed);
        try self.append(rec.qual);
        try self.append(rec.aux_raw);

        try self.flushBlocks(false);
    }

    /// Flush the final partial block and the BGZF EOF marker.
    pub fn finish(self: *Writer) !void {
        try self.flushBlocks(true);
        try writeEof(self.out);
    }

    pub fn deinit(self: *Writer) void {
        self.buf.deinit(self.gpa);
        self.* = undefined;
    }

    fn append(self: *Writer, bytes: []const u8) !void {
        try self.buf.appendSlice(self.gpa, bytes);
    }
    fn appendByte(self: *Writer, b: u8) !void {
        try self.buf.append(self.gpa, b);
    }
    fn appendU16(self: *Writer, v: u16) !void {
        var t: [2]u8 = undefined;
        std.mem.writeInt(u16, &t, v, .little);
        try self.append(&t);
    }
    fn appendU32(self: *Writer, v: u32) !void {
        var t: [4]u8 = undefined;
        std.mem.writeInt(u32, &t, v, .little);
        try self.append(&t);
    }
    fn appendI32(self: *Writer, v: i32) !void {
        var t: [4]u8 = undefined;
        std.mem.writeInt(i32, &t, v, .little);
        try self.append(&t);
    }

    fn flushBlocks(self: *Writer, final: bool) !void {
        while (self.buf.items.len >= max_block_payload) {
            try writeBlock(self.out, self.buf.items[0..max_block_payload]);
            const rem = self.buf.items.len - max_block_payload;
            std.mem.copyForwards(u8, self.buf.items[0..rem], self.buf.items[max_block_payload..]);
            self.buf.items.len = rem;
        }
        if (final) {
            if (self.buf.items.len > 0) try writeBlock(self.out, self.buf.items);
            self.buf.items.len = 0;
        }
    }
};

pub fn writer(gpa: Allocator, out: *std.Io.Writer) Writer {
    return .{ .gpa = gpa, .out = out, .buf = .empty };
}

// Tests ---------------------------------------------------------------------

test "bgzf writeBlock/inflateAll round-trip" {
    const payload = "BGZF block round-trip payload: ACGTACGT the quick brown fox 0123456789, repeated content compresses well ACGTACGT ACGTACGT.";
    var out: std.Io.Writer.Allocating = .init(testing.allocator);
    defer out.deinit();
    try writeBlock(&out.writer, payload);
    try writeEof(&out.writer);

    var in: std.Io.Reader = .fixed(out.written());
    const back = try inflateAll(testing.allocator, &in);
    defer testing.allocator.free(back);
    try testing.expectEqualStrings(payload, back);
}

test "read embedded t1.bam" {
    const file = @embedFile("./test/t1.bam");
    var in: std.Io.Reader = .fixed(file);
    var br = try open(testing.allocator, &in);
    defer br.deinit();

    try testing.expectEqual(@as(usize, 4), br.header.refs.len);
    try testing.expectEqual(@as(?i32, 0), br.header.tid("insert"));
    try testing.expectEqual(@as(?i32, 1), br.header.tid("ref1"));
    try testing.expectEqual(@as(?u64, 45), br.header.refLen(1));

    var count: usize = 0;
    while (try br.next()) |_| count += 1;
    try testing.expectEqual(@as(usize, 15), count);
}

test "SAM/BAM cross-check on identical content" {
    const sam_text = @embedFile("./test/t1_bam.sam");
    var sam_in: std.Io.Reader = .fixed(sam_text);
    var sr = sam.reader(testing.allocator, &sam_in);
    var sh = try sr.readHeader();
    defer sh.deinit(testing.allocator);

    const bam_bytes = @embedFile("./test/t1.bam");
    var bam_in: std.Io.Reader = .fixed(bam_bytes);
    var br = try open(testing.allocator, &bam_in);
    defer br.deinit();

    while (try sr.next(&sh)) |sa| {
        const ba = (try br.next()) orelse return error.TestUnexpectedResult;

        try testing.expectEqualStrings(sa.qname, ba.qname);
        try testing.expectEqual(sa.flag, ba.flag);
        try testing.expectEqual(@as(i64, sa.pos), ba.pos);
        try testing.expectEqual(sa.mapq, ba.mapq);

        if (sa.tid >= 0) {
            try testing.expectEqualStrings(sh.refName(sa.tid).?, br.header.refName(ba.ref_id).?);
        } else {
            try testing.expectEqual(@as(i32, -1), ba.ref_id);
        }

        var sc = sa.cigar();
        var bc = ba.cigar();
        while (true) {
            const s = try sc.next();
            const b = bc.next();
            if (s == null and b == null) break;
            try testing.expect(s != null and b != null);
            try testing.expectEqual(s.?.op, b.?.op);
            try testing.expectEqual(s.?.len, b.?.len);
        }

        if (sa.seq.len > 0) {
            try testing.expectEqual(@as(u32, @intCast(sa.seq.len)), ba.l_seq);
            for (sa.seq, 0..) |ch, i| try testing.expectEqual(ch, ba.seqBase(i));
        }
        if (sa.qual.len > 0) {
            try testing.expectEqual(sa.qual.len, ba.qual.len);
            for (sa.qual, 0..) |q, i| try testing.expectEqual(q, ba.qual[i] + 33);
        }
    }
    try testing.expect((try br.next()) == null);
}

test "BAM round-trip: read -> write -> read" {
    const bam_bytes = @embedFile("./test/t1.bam");
    var in1: std.Io.Reader = .fixed(bam_bytes);
    var br1 = try open(testing.allocator, &in1);
    defer br1.deinit();

    var recs: std.ArrayList(Record) = .empty;
    defer recs.deinit(testing.allocator);
    while (try br1.next()) |r| try recs.append(testing.allocator, r);

    var out: std.Io.Writer.Allocating = .init(testing.allocator);
    defer out.deinit();
    var bw = writer(testing.allocator, &out.writer);
    defer bw.deinit();
    try bw.writeHeader(br1.header_text, br1.header.refs);
    for (recs.items) |*r| try bw.writeRecord(r);
    try bw.finish();

    var in2: std.Io.Reader = .fixed(out.written());
    var br2 = try open(testing.allocator, &in2);
    defer br2.deinit();

    try testing.expectEqual(br1.header.refs.len, br2.header.refs.len);
    var i: usize = 0;
    while (try br2.next()) |r2| : (i += 1) {
        const r1 = recs.items[i];
        try testing.expectEqualStrings(r1.qname, r2.qname);
        try testing.expectEqual(r1.flag, r2.flag);
        try testing.expectEqual(r1.ref_id, r2.ref_id);
        try testing.expectEqual(r1.pos, r2.pos);
        try testing.expectEqual(r1.mapq, r2.mapq);
        try testing.expectEqual(r1.bin, r2.bin);
        try testing.expectEqual(r1.next_ref_id, r2.next_ref_id);
        try testing.expectEqual(r1.next_pos, r2.next_pos);
        try testing.expectEqual(r1.tlen, r2.tlen);
        try testing.expectEqual(r1.l_seq, r2.l_seq);
        try testing.expectEqualSlices(u8, r1.cigar_raw, r2.cigar_raw);
        try testing.expectEqualSlices(u8, r1.seq_packed, r2.seq_packed);
        try testing.expectEqualSlices(u8, r1.qual, r2.qual);
        try testing.expectEqualSlices(u8, r1.aux_raw, r2.aux_raw);
    }
    try testing.expectEqual(recs.items.len, i);
}

test "seqNt16 encode/decode round-trip" {
    const bases = "=ACMGRSVTWYHKDBN";
    for (bases, 0..) |ch, code| {
        try testing.expectEqual(@as(u4, @intCast(code)), seqNt16(ch));
        try testing.expectEqual(ch, seq_nt16_str[code]);
    }
}
