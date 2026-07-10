const std = @import("std");
const testing = std.testing;

// Sequence alphabets, modeled on biogo's `alphabet` package. Provides:
//   - byte-based `Alphabet` values (DNA/RNA/protein, + gapped/IUPAC-redundant)
//   - complement / reverse-complement for the nucleic alphabets
//   - per-letter enums `Nucleotide` (IUPAC 4-bit) and `AminoAcid`
//   - FASTQ quality `Encoding` (Sanger/Solexa/Illumina) + Phred/Solexa helpers

pub const MolType = enum { dna, rna, protein };

// Byte-based alphabets ------------------------------------------------------

pub const Alphabet = struct {
    letters: []const u8,
    valid: [256]bool,
    index: [256]i16, // letter -> position in `letters`, -1 if invalid
    gap: u8,
    ambiguous: u8,
    case_sensitive: bool,
    mol_type: MolType,
    complement_table: ?*const [256]u8, // null for non-nucleic

    pub fn len(self: *const Alphabet) usize {
        return self.letters.len;
    }
    pub fn isValid(self: *const Alphabet, c: u8) bool {
        return self.valid[c];
    }
    pub fn indexOf(self: *const Alphabet, c: u8) i16 {
        return self.index[c];
    }
    pub fn letter(self: *const Alphabet, i: usize) u8 {
        return self.letters[i];
    }
    pub fn isGap(self: *const Alphabet, c: u8) bool {
        return c == self.gap;
    }
    pub fn isCased(self: *const Alphabet) bool {
        return self.case_sensitive;
    }
    /// Index of the first invalid letter in `seq`, or null if all are valid.
    pub fn allValid(self: *const Alphabet, seq: []const u8) ?usize {
        for (seq, 0..) |c, i| {
            if (!self.valid[c]) return i;
        }
        return null;
    }
    /// Complement of a single letter, or null if this alphabet has no
    /// complement (protein) or the letter has none.
    pub fn complement(self: *const Alphabet, c: u8) ?u8 {
        const table = self.complement_table orelse return null;
        const r = table[c];
        return if (r == 0) null else r;
    }
    /// Write the complement of each `src` byte into `dest` (same order). Bytes
    /// without a complement pass through unchanged. `dest.len` must be >= src.len.
    pub fn complementInto(self: *const Alphabet, dest: []u8, src: []const u8) void {
        for (src, 0..) |c, i| dest[i] = self.complement(c) orelse c;
    }
    /// Write the reverse-complement of `src` into `dest`.
    pub fn reverseComplementInto(self: *const Alphabet, dest: []u8, src: []const u8) void {
        const n = src.len;
        for (src, 0..) |c, i| dest[n - 1 - i] = self.complement(c) orelse c;
    }
    /// Reverse-complement `seq` in place.
    pub fn reverseComplementInPlace(self: *const Alphabet, seq: []u8) void {
        var i: usize = 0;
        var j: usize = seq.len;
        while (i < j) {
            j -= 1;
            const a = self.complement(seq[i]) orelse seq[i];
            const b = self.complement(seq[j]) orelse seq[j];
            seq[i] = b;
            seq[j] = a;
            i += 1;
        }
    }
};

fn otherCase(c: u8) ?u8 {
    return switch (c) {
        'a'...'z' => c - 32,
        'A'...'Z' => c + 32,
        else => null,
    };
}

fn buildAlphabet(
    comptime letters: []const u8,
    comptime mol: MolType,
    comptime gap: u8,
    comptime ambiguous: u8,
    comptime case_sensitive: bool,
    comptime comp: ?*const [256]u8,
) Alphabet {
    comptime {
        var valid: [256]bool = undefined;
        var index: [256]i16 = undefined;
        for (&valid) |*v| v.* = false;
        for (&index) |*v| v.* = -1;
        for (letters, 0..) |c, i| {
            valid[c] = true;
            index[c] = @intCast(i);
            if (!case_sensitive) {
                if (otherCase(c)) |oc| {
                    valid[oc] = true;
                    index[oc] = @intCast(i);
                }
            }
        }
        return .{
            .letters = letters,
            .valid = valid,
            .index = index,
            .gap = gap,
            .ambiguous = ambiguous,
            .case_sensitive = case_sensitive,
            .mol_type = mol,
            .complement_table = comp,
        };
    }
}

fn buildComplement(comptime from: []const u8, comptime to: []const u8) [256]u8 {
    comptime {
        var t: [256]u8 = undefined;
        for (&t) |*v| v.* = 0; // 0 = no complement
        for (from, to) |f, c| t[f] = c;
        return t;
    }
}

// biogo complement pairings (verbatim letter sets).
const dna_comp = buildComplement("acgtnxACGTNX-", "tgcanxTGCANX-");
const dna_red_comp = buildComplement(
    "acmgrsvtwyhkdbnxACMGRSVTWYHKDBNX-",
    "tgkcysbawrdmhvnxTGKCYSBAWRDMHVNX-",
);
const rna_comp = buildComplement("acgunxACGUNX-", "ugcanxUGCANX-");
const rna_red_comp = buildComplement(
    "acmgrsvuwyhkdbnxACMGRSVUWYHKDBNX-",
    "ugkcysbawrdmhvnxUGKCYSBAWRDMHVNX-",
);

pub const dna = buildAlphabet("acgt", .dna, '-', 'n', false, &dna_comp);
pub const dna_gapped = buildAlphabet("-acgt", .dna, '-', 'n', false, &dna_comp);
pub const dna_redundant = buildAlphabet("-acmgrsvtwyhkdbn", .dna, '-', 'n', false, &dna_red_comp);
pub const rna = buildAlphabet("acgu", .rna, '-', 'n', false, &rna_comp);
pub const rna_gapped = buildAlphabet("-acgu", .rna, '-', 'n', false, &rna_comp);
pub const rna_redundant = buildAlphabet("-acmgrsvuwyhkdbn", .rna, '-', 'n', false, &rna_red_comp);
pub const protein = buildAlphabet("-abcdefghijklmnpqrstvwxyz*", .protein, '-', 'x', false, null);

pub const Kind = enum {
    dna,
    dna_gapped,
    dna_redundant,
    rna,
    rna_gapped,
    rna_redundant,
    protein,

    pub fn alphabet(self: Kind) *const Alphabet {
        return switch (self) {
            .dna => &dna,
            .dna_gapped => &dna_gapped,
            .dna_redundant => &dna_redundant,
            .rna => &rna,
            .rna_gapped => &rna_gapped,
            .rna_redundant => &rna_redundant,
            .protein => &protein,
        };
    }
    pub fn molType(self: Kind) MolType {
        return self.alphabet().mol_type;
    }
};

// Per-letter enums ----------------------------------------------------------

/// IUPAC nucleotide, encoded as the 4-bit code used by BAM (`=ACMGRSVTWYHKDBN`).
pub const Nucleotide = enum(u4) {
    eq = 0, // '='
    a = 1,
    c = 2,
    m = 3,
    g = 4,
    r = 5,
    s = 6,
    v = 7,
    t = 8,
    w = 9,
    y = 10,
    h = 11,
    k = 12,
    d = 13,
    b = 14,
    n = 15,

    pub const chars = "=ACMGRSVTWYHKDBN";

    pub fn fromChar(ch: u8) ?Nucleotide {
        return switch (ch) {
            '=' => .eq,
            'A', 'a' => .a,
            'C', 'c' => .c,
            'M', 'm' => .m,
            'G', 'g' => .g,
            'R', 'r' => .r,
            'S', 's' => .s,
            'V', 'v' => .v,
            'T', 't' => .t,
            'U', 'u' => .t, // uracil shares the T code (BAM convention)
            'W', 'w' => .w,
            'Y', 'y' => .y,
            'H', 'h' => .h,
            'K', 'k' => .k,
            'D', 'd' => .d,
            'B', 'b' => .b,
            'N', 'n' => .n,
            else => null,
        };
    }

    pub fn toChar(self: Nucleotide) u8 {
        return chars[@intFromEnum(self)];
    }

    /// IUPAC complement within the 16-code space.
    pub fn complement(self: Nucleotide) Nucleotide {
        return switch (self) {
            .eq => .eq,
            .a => .t,
            .t => .a,
            .c => .g,
            .g => .c,
            .m => .k,
            .k => .m,
            .r => .y,
            .y => .r,
            .s => .s,
            .w => .w,
            .n => .n,
            .v => .b,
            .b => .v,
            .h => .d,
            .d => .h,
        };
    }
};

/// Amino acid: 20 standard + Sec/Pyl + ambiguity (Asx/Glx/Xle) + X/stop/gap.
pub const AminoAcid = enum(u8) {
    ala, // A
    arg, // R
    asn, // N
    asp, // D
    cys, // C
    gln, // Q
    glu, // E
    gly, // G
    his, // H
    ile, // I
    leu, // L
    lys, // K
    met, // M
    phe, // F
    pro, // P
    ser, // S
    thr, // T
    trp, // W
    tyr, // Y
    val, // V
    sec, // U selenocysteine
    pyl, // O pyrrolysine
    asx, // B Asn or Asp
    glx, // Z Gln or Glu
    xle, // J Leu or Ile
    unknown, // X
    stop, // *
    gap, // -

    pub const one_letter = "ARNDCQEGHILKMFPSTWYVUOBZJX*-";

    pub fn toChar(self: AminoAcid) u8 {
        return one_letter[@intFromEnum(self)];
    }

    pub fn fromChar(ch: u8) ?AminoAcid {
        const u = if (ch >= 'a' and ch <= 'z') ch - 32 else ch;
        return switch (u) {
            'A' => .ala,
            'R' => .arg,
            'N' => .asn,
            'D' => .asp,
            'C' => .cys,
            'Q' => .gln,
            'E' => .glu,
            'G' => .gly,
            'H' => .his,
            'I' => .ile,
            'L' => .leu,
            'K' => .lys,
            'M' => .met,
            'F' => .phe,
            'P' => .pro,
            'S' => .ser,
            'T' => .thr,
            'W' => .trp,
            'Y' => .tyr,
            'V' => .val,
            'U' => .sec,
            'O' => .pyl,
            'B' => .asx,
            'Z' => .glx,
            'J' => .xle,
            'X' => .unknown,
            '*' => .stop,
            '-' => .gap,
            else => null,
        };
    }

    pub fn threeLetter(self: AminoAcid) []const u8 {
        return switch (self) {
            .ala => "Ala",
            .arg => "Arg",
            .asn => "Asn",
            .asp => "Asp",
            .cys => "Cys",
            .gln => "Gln",
            .glu => "Glu",
            .gly => "Gly",
            .his => "His",
            .ile => "Ile",
            .leu => "Leu",
            .lys => "Lys",
            .met => "Met",
            .phe => "Phe",
            .pro => "Pro",
            .ser => "Ser",
            .thr => "Thr",
            .trp => "Trp",
            .tyr => "Tyr",
            .val => "Val",
            .sec => "Sec",
            .pyl => "Pyl",
            .asx => "Asx",
            .glx => "Glx",
            .xle => "Xle",
            .unknown => "Xaa",
            .stop => "Ter",
            .gap => "Gap",
        };
    }
};

// Quality encodings ---------------------------------------------------------

pub const Phred = u8; // Phred quality score
pub const Solexa = i8; // Solexa quality score

/// Error probability for a Phred score: 10^(-q/10).
pub fn phredErrorProb(q: Phred) f64 {
    return std.math.pow(f64, 10.0, -@as(f64, @floatFromInt(q)) / 10.0);
}

/// Nearest Phred score for an error probability.
pub fn phredFromErrorProb(p: f64) Phred {
    return clampPhred(-10.0 * std.math.log10(p));
}

/// Solexa -> Phred: Q_phred = 10*log10(10^(Q_sol/10) + 1).
pub fn solexaToPhred(s: Solexa) Phred {
    const p = std.math.pow(f64, 10.0, @as(f64, @floatFromInt(s)) / 10.0);
    return clampPhred(10.0 * std.math.log10(p + 1.0));
}

/// Phred -> Solexa: Q_sol = 10*log10(10^(Q_phred/10) - 1).
pub fn phredToSolexa(q: Phred) Solexa {
    if (q == 0) return -5; // 10^0 - 1 = 0 -> -inf; biogo's practical floor
    const p = std.math.pow(f64, 10.0, @as(f64, @floatFromInt(q)) / 10.0);
    return clampSolexa(10.0 * std.math.log10(p - 1.0));
}

fn clampPhred(x: f64) Phred {
    const r = @round(x);
    if (r <= 0) return 0;
    if (r >= 255) return 255;
    return @intFromFloat(r);
}

fn clampSolexa(x: f64) Solexa {
    const r = @round(x);
    if (r <= -128) return -128;
    if (r >= 127) return 127;
    return @intFromFloat(r);
}

/// FASTQ quality encoding schemes.
pub const Encoding = enum {
    sanger, // Phred+33
    solexa, // Solexa+64
    illumina_1_3, // Phred+64
    illumina_1_5, // Phred+64
    illumina_1_8, // Phred+33

    /// ASCII offset of the scheme.
    pub fn offset(self: Encoding) u8 {
        return switch (self) {
            .sanger, .illumina_1_8 => 33,
            .solexa, .illumina_1_3, .illumina_1_5 => 64,
        };
    }

    /// Decode an ASCII quality byte to a Phred score.
    pub fn decodeToPhred(self: Encoding, ascii: u8) Phred {
        return switch (self) {
            .sanger, .illumina_1_8 => ascii - 33,
            .illumina_1_3, .illumina_1_5 => ascii - 64,
            .solexa => solexaToPhred(@intCast(@as(i16, ascii) - 64)),
        };
    }

    /// Encode a Phred score to an ASCII quality byte.
    pub fn encodePhred(self: Encoding, q: Phred) u8 {
        return switch (self) {
            .sanger, .illumina_1_8 => q + 33,
            .illumina_1_3, .illumina_1_5 => q + 64,
            .solexa => @intCast(@as(i16, phredToSolexa(q)) + 64),
        };
    }
};

// Tests ---------------------------------------------------------------------

test "alphabet: validity, index, gap, ambiguous" {
    try testing.expect(dna.isValid('A'));
    try testing.expect(dna.isValid('a'));
    try testing.expect(dna.isValid('t'));
    try testing.expect(!dna.isValid('x'));
    try testing.expect(!dna.isValid('n')); // ambiguous letter not in plain DNA set
    try testing.expect(dna_redundant.isValid('n'));

    try testing.expectEqual(@as(usize, 4), dna.len());
    try testing.expectEqual(@as(i16, 0), dna.indexOf('a'));
    try testing.expectEqual(@as(i16, 0), dna.indexOf('A'));
    try testing.expectEqual(@as(u8, 'a'), dna.letter(0));
    try testing.expectEqual(@as(i16, -1), dna.indexOf('x'));

    try testing.expectEqual(@as(u8, '-'), dna.gap);
    try testing.expectEqual(@as(u8, 'n'), dna.ambiguous);
    try testing.expect(dna_gapped.isGap('-'));

    try testing.expectEqual(@as(?usize, null), dna.allValid("ACGTacgt"));
    try testing.expectEqual(@as(?usize, 4), dna.allValid("ACGTxA"));
}

test "Kind dispatch" {
    try testing.expectEqual(MolType.dna, Kind.dna.molType());
    try testing.expectEqual(MolType.rna, Kind.rna.molType());
    try testing.expectEqual(MolType.protein, Kind.protein.molType());
    try testing.expect(Kind.rna.alphabet().isValid('u'));
    try testing.expect(!Kind.dna.alphabet().isValid('u'));
}

test "complement and reverse-complement" {
    try testing.expectEqual(@as(?u8, 't'), dna.complement('a'));
    try testing.expectEqual(@as(?u8, 'A'), dna.complement('T'));
    try testing.expectEqual(@as(?u8, null), protein.complement('a'));
    try testing.expectEqual(@as(?u8, 'y'), dna_redundant.complement('r')); // IUPAC

    var buf: [4]u8 = undefined;
    dna.reverseComplementInto(&buf, "AACG");
    try testing.expectEqualStrings("CGTT", &buf);

    var seq = [_]u8{ 'A', 'A', 'C', 'G' };
    dna.reverseComplementInPlace(&seq);
    try testing.expectEqualStrings("CGTT", &seq);

    var odd = [_]u8{ 'A', 'C', 'G' };
    dna.reverseComplementInPlace(&odd);
    try testing.expectEqualStrings("CGT", &odd);
}

test "Nucleotide enum matches seq_nt16 ordering and complements" {
    const chars = "=ACMGRSVTWYHKDBN";
    for (chars, 0..) |ch, code| {
        const nuc = Nucleotide.fromChar(ch).?;
        try testing.expectEqual(@as(u4, @intCast(code)), @intFromEnum(nuc));
        try testing.expectEqual(ch, nuc.toChar());
        if (ch >= 'A' and ch <= 'Z') {
            try testing.expectEqual(nuc, Nucleotide.fromChar(ch + 32).?);
        }
    }
    try testing.expectEqual(Nucleotide.t, Nucleotide.fromChar('U').?);
    try testing.expectEqual(Nucleotide.t, Nucleotide.a.complement());
    try testing.expectEqual(Nucleotide.a, Nucleotide.t.complement());
    try testing.expectEqual(Nucleotide.y, Nucleotide.r.complement());
    try testing.expectEqual(@as(?Nucleotide, null), Nucleotide.fromChar('z'));
}

test "AminoAcid enum" {
    try testing.expectEqual(AminoAcid.trp, AminoAcid.fromChar('W').?);
    try testing.expectEqual(AminoAcid.trp, AminoAcid.fromChar('w').?);
    try testing.expectEqual(@as(u8, 'W'), AminoAcid.trp.toChar());
    try testing.expectEqualStrings("Trp", AminoAcid.trp.threeLetter());
    try testing.expectEqual(AminoAcid.stop, AminoAcid.fromChar('*').?);
    try testing.expectEqualStrings("Ter", AminoAcid.stop.threeLetter());
    try testing.expectEqual(@as(?AminoAcid, null), AminoAcid.fromChar('!'));

    for (AminoAcid.one_letter, 0..) |ch, code| {
        const aa: AminoAcid = @enumFromInt(code);
        try testing.expectEqual(ch, aa.toChar());
        try testing.expectEqual(aa, AminoAcid.fromChar(ch).?);
    }
}

test "quality encodings" {
    try testing.expectEqual(@as(Phred, 0), Encoding.sanger.decodeToPhred('!'));
    try testing.expectEqual(@as(Phred, 40), Encoding.sanger.decodeToPhred('I'));
    try testing.expectEqual(@as(u8, '!'), Encoding.sanger.encodePhred(0));
    try testing.expectEqual(@as(u8, 'I'), Encoding.sanger.encodePhred(40));

    try testing.expectEqual(@as(u8, 64), Encoding.illumina_1_3.offset());
    try testing.expectEqual(@as(Phred, 30), Encoding.illumina_1_3.decodeToPhred(30 + 64));

    try testing.expectApproxEqAbs(@as(f64, 0.01), phredErrorProb(20), 1e-9);
    try testing.expectApproxEqAbs(@as(f64, 0.001), phredErrorProb(30), 1e-9);
    try testing.expectEqual(@as(Phred, 20), phredFromErrorProb(0.01));

    // High Phred scores round-trip through Solexa.
    try testing.expectEqual(@as(Phred, 30), solexaToPhred(phredToSolexa(30)));
}

// Rank transform & q-grams --------------------------------------------------

/// Maps alphabet symbols to their lexicographic rank (ascending byte order) and
/// enumerates rank-packed q-grams. Ported from rust-bio's `alphabets::RankTransform`.
/// This is distinct from `Alphabet.index` (position within a fixed `letters`
/// string): ranks here are assigned by ascending byte value over the symbol set.
pub const RankTransform = struct {
    ranks: [256]u8 = @splat(0),
    present: [256]bool = @splat(false),
    len: usize = 0,

    /// Build a transform over the given symbol bytes. Ranks are assigned in
    /// ascending byte order (0, 1, 2, ...).
    pub fn init(symbols: []const u8) RankTransform {
        var rt = RankTransform{};
        for (symbols) |s| rt.present[s] = true;
        var r: u8 = 0;
        var c: usize = 0;
        while (c < 256) : (c += 1) {
            if (rt.present[c]) {
                rt.ranks[c] = r;
                r += 1;
            }
        }
        rt.len = r;
        return rt;
    }

    /// Rank of symbol `a`. Asserts `a` is in the alphabet.
    pub fn get(self: *const RankTransform, a: u8) u8 {
        std.debug.assert(self.present[a]);
        return self.ranks[a];
    }

    /// Bits needed to encode the largest rank: ceil(log2(len)). Suitable as the
    /// `width` for `bitenc.BitEnc`.
    pub fn getWidth(self: *const RankTransform) u6 {
        if (self.len <= 1) return 0;
        return @intCast(std.math.log2_int_ceil(usize, self.len));
    }

    /// Transform `text` into an owned slice of ranks. Caller frees.
    pub fn transform(self: *const RankTransform, gpa: std.mem.Allocator, text: []const u8) ![]u8 {
        const out = try gpa.alloc(u8, text.len);
        for (text, 0..) |c, i| out[i] = self.get(c);
        return out;
    }

    /// Iterate q-grams of `text`, each encoded as a usize by packing symbol
    /// ranks in `getWidth()` bits. Asserts `q > 0` and `q * width <= 64`.
    pub fn qgrams(self: *const RankTransform, q: u32, text: []const u8) QGrams {
        std.debug.assert(q > 0);
        const bits = self.getWidth();
        const shift: usize = @as(usize, bits) * q;
        std.debug.assert(shift <= @bitSizeOf(usize));
        const mask: usize = if (shift >= @bitSizeOf(usize))
            std.math.maxInt(usize)
        else
            (@as(usize, 1) << @intCast(shift)) - 1;
        var g = QGrams{ .ranks = self, .text = text, .pos = 0, .bits = bits, .mask = mask, .qgram = 0 };
        var i: u32 = 0;
        while (i < q - 1) : (i += 1) _ = g.next();
        return g;
    }
};

/// Iterator over rank-packed q-grams. See `RankTransform.qgrams`.
pub const QGrams = struct {
    ranks: *const RankTransform,
    text: []const u8,
    pos: usize,
    bits: u6,
    mask: usize,
    qgram: usize,

    pub fn next(self: *QGrams) ?usize {
        if (self.pos >= self.text.len) return null;
        const c = self.text[self.pos];
        self.pos += 1;
        const b = self.ranks.get(c);
        self.qgram = ((self.qgram << self.bits) | b) & self.mask;
        return self.qgram;
    }
};

test "RankTransform transform" {
    const rt = RankTransform.init("ACGTacgt");
    const out = try rt.transform(testing.allocator, "aAcCgGtT");
    defer testing.allocator.free(out);
    try testing.expectEqualSlices(u8, &.{ 4, 0, 5, 1, 6, 2, 7, 3 }, out);
    try testing.expectEqual(@as(u8, 0), rt.get('A'));
    try testing.expectEqual(@as(u8, 7), rt.get('t'));
}

test "RankTransform getWidth" {
    try testing.expectEqual(@as(u6, 2), RankTransform.init("ACGT").getWidth());
    try testing.expectEqual(@as(u6, 3), RankTransform.init("ACGTN").getWidth());
}

test "RankTransform qgrams" {
    const rt = RankTransform.init("ACGTacgt");
    var g = rt.qgrams(2, "ACGT");
    try testing.expectEqual(@as(?usize, 1), g.next());
    try testing.expectEqual(@as(?usize, 10), g.next());
    try testing.expectEqual(@as(?usize, 19), g.next());
    try testing.expectEqual(@as(?usize, null), g.next());
}

test "RankTransform qgram shift-left no overflow (q*bits == 64)" {
    const rt = RankTransform.init("ACTG");
    var buf: [400]u8 = undefined;
    for (0..100) |i| @memcpy(buf[i * 4 ..][0..4], "ACTG");
    var g = rt.qgrams(@bitSizeOf(usize) / 2, &buf); // q=32, bits=2
    _ = g.next(); // must not panic on the 1<<64 mask edge
}
