//! Benchmarks for algorithms ported from rust-bio. Grown one tranche at a time.
//!
//! Fixed inputs are built once in `register()` (which runs before timing) and
//! stored at module scope — zBench times the whole body of a BenchFunc, so keep
//! per-run setup out of it.

const std = @import("std");
const zbench = @import("zbench");
const zbio = @import("zbio");

const N = 100_000;
var seq_a: [N]u8 = undefined;
var seq_b: [N]u8 = undefined;

// A rank transform over DNA, built once for the qgram bench.
var dna_ranks: zbio.alphabet.RankTransform = undefined;

pub fn register(bench: *zbench.Benchmark) !void {
    fillSeq(&seq_a, 1);
    fillSeq(&seq_b, 2);
    dna_ranks = zbio.alphabet.RankTransform.init("ACGT");

    // Tranche 1 (foundational).
    try bench.add("hamming 100k", benchHamming, .{ .items_per_run = N });
    try bench.add("levenshtein 512x512", benchLevenshtein, .{});
    try bench.add("gc content 100k", benchGcContent, .{ .bytes_per_run = N });
    try bench.add("qgrams(8) 100k", benchQgrams, .{ .bytes_per_run = N });

    // Tranche 2 (pairwise alignment).
    try bench.add("global SW 512x512", benchGlobalSW, .{});
    try bench.add("local SW 512x512", benchLocalSW, .{});
    try bench.add("local blosum62 protein", benchBlosumLocal, .{});

    // Tranche 3 (exact pattern matching): scan a 12-mer over the 100k haystack.
    try bench.add("horspool search 100k", benchHorspool, .{ .bytes_per_run = N });
    try bench.add("shift-and search 100k", benchShiftAnd, .{ .bytes_per_run = N });
    try bench.add("bndm search 100k", benchBndm, .{ .bytes_per_run = N });
    try bench.add("kmp search 100k", benchKmp, .{ .bytes_per_run = N });
    try bench.add("bom search 100k", benchBom, .{ .bytes_per_run = N });

    // Tranche 4 (FM-index stack). SA build is O(n log^2 n) here, so use 4k bases.
    sa_text[SA_N] = '$';
    @memcpy(sa_text[0..SA_N], seq_a[0..SA_N]);
    try bench.add("suffix array build 4k", benchSaBuild, .{ .bytes_per_run = SA_N });
    try bench.add("bwt build 4k", benchBwtBuild, .{ .bytes_per_run = SA_N });
    try bench.add("fm-index build 4k", benchFmBuild, .{ .track_allocations = true });
    try bench.add("fm-index search", benchFmSearch, .{});
}

const pm = zbio.pm;
const PAT = "ACGTACGTACGT";

fn countMatches(it_ptr: anytype) usize {
    var n: usize = 0;
    while (it_ptr.next()) |_| n +%= 1;
    return n;
}

fn benchHorspool(_: std.mem.Allocator) void {
    const hs = pm.horspool.Horspool.init(PAT);
    var it = hs.findAll(&seq_a);
    std.mem.doNotOptimizeAway(countMatches(&it));
}

fn benchShiftAnd(_: std.mem.Allocator) void {
    const sa = pm.shift_and.ShiftAnd.init(PAT);
    var it = sa.findAll(&seq_a);
    std.mem.doNotOptimizeAway(countMatches(&it));
}

fn benchBndm(_: std.mem.Allocator) void {
    const b = pm.bndm.BNDM.init(PAT);
    var it = b.findAll(&seq_a);
    std.mem.doNotOptimizeAway(countMatches(&it));
}

fn benchKmp(gpa: std.mem.Allocator) void {
    var k = pm.kmp.KMP.init(gpa, PAT) catch |e| std.debug.panic("kmp: {t}", .{e});
    defer k.deinit(gpa);
    var it = k.findAll(&seq_a);
    std.mem.doNotOptimizeAway(countMatches(&it));
}

fn benchBom(gpa: std.mem.Allocator) void {
    var b = pm.bom.BOM.init(gpa, PAT) catch |e| std.debug.panic("bom: {t}", .{e});
    defer b.deinit(gpa);
    var it = b.findAll(&seq_a);
    std.mem.doNotOptimizeAway(countMatches(&it));
}

const idx = zbio.index;
const SA_N = 4096;
var sa_text: [SA_N + 1]u8 = undefined; // seq_a[0..SA_N] ++ '$'

fn benchSaBuild(gpa: std.mem.Allocator) void {
    const sa = idx.suffix_array.suffixArray(gpa, &sa_text) catch |e| std.debug.panic("sa: {t}", .{e});
    defer gpa.free(sa);
    std.mem.doNotOptimizeAway(sa[0]);
}

fn benchBwtBuild(gpa: std.mem.Allocator) void {
    const sa = idx.suffix_array.suffixArray(gpa, &sa_text) catch |e| std.debug.panic("sa: {t}", .{e});
    defer gpa.free(sa);
    const b = idx.bwt.bwt(gpa, &sa_text, sa) catch |e| std.debug.panic("bwt: {t}", .{e});
    defer gpa.free(b);
    std.mem.doNotOptimizeAway(b[0]);
}

fn benchFmBuild(gpa: std.mem.Allocator) void {
    const sa = idx.suffix_array.suffixArray(gpa, &sa_text) catch |e| std.debug.panic("sa: {t}", .{e});
    defer gpa.free(sa);
    const b = idx.bwt.bwt(gpa, &sa_text, sa) catch |e| std.debug.panic("bwt: {t}", .{e});
    defer gpa.free(b);
    const l = idx.bwt.less(gpa, b) catch |e| std.debug.panic("less: {t}", .{e});
    defer gpa.free(l);
    var occ = idx.bwt.Occ.init(gpa, b, 32) catch |e| std.debug.panic("occ: {t}", .{e});
    defer occ.deinit(gpa);
    std.mem.doNotOptimizeAway(l[0]);
}

// Lazily build a persistent FM-index (leaked into the process) so the search
// benchmark measures pure query time, not construction.
var fm_ready = false;
var fm_bwt: []u8 = undefined;
var fm_less: []usize = undefined;
var fm_occ: idx.bwt.Occ = undefined;
var fm_index: idx.fmindex.FMIndex = undefined;

fn ensureFm() void {
    if (fm_ready) return;
    const a = std.heap.page_allocator;
    const sa = idx.suffix_array.suffixArray(a, &sa_text) catch |e| std.debug.panic("sa: {t}", .{e});
    fm_bwt = idx.bwt.bwt(a, &sa_text, sa) catch |e| std.debug.panic("bwt: {t}", .{e});
    a.free(sa);
    fm_less = idx.bwt.less(a, fm_bwt) catch |e| std.debug.panic("less: {t}", .{e});
    fm_occ = idx.bwt.Occ.init(a, fm_bwt, 32) catch |e| std.debug.panic("occ: {t}", .{e});
    fm_index = idx.fmindex.FMIndex.init(fm_bwt, fm_less, &fm_occ);
    fm_ready = true;
}

fn benchFmSearch(_: std.mem.Allocator) void {
    ensureFm();
    const bsr = fm_index.backwardSearch("ACGTACGT");
    const found: usize = switch (bsr) {
        .complete => |iv| iv.upper - iv.lower,
        .partial => |p| p.interval.upper - p.interval.lower,
        .absent => 0,
    };
    std.mem.doNotOptimizeAway(found);
}

const pw = zbio.alignment.pairwise;

fn benchGlobalSW(gpa: std.mem.Allocator) void {
    var a = pw.Aligner(pw.MatchParams).init(gpa, -5, -1, pw.MatchParams.init(1, -1));
    defer a.deinit();
    var aln = a.global(seq_a[0..512], seq_b[0..512]) catch |e| std.debug.panic("sw: {t}", .{e});
    defer aln.deinit(gpa);
    std.mem.doNotOptimizeAway(aln.score);
}

fn benchLocalSW(gpa: std.mem.Allocator) void {
    var a = pw.Aligner(pw.MatchParams).init(gpa, -5, -1, pw.MatchParams.init(1, -1));
    defer a.deinit();
    var aln = a.local(seq_a[0..512], seq_b[0..512]) catch |e| std.debug.panic("sw: {t}", .{e});
    defer aln.deinit(gpa);
    std.mem.doNotOptimizeAway(aln.score);
}

fn benchBlosumLocal(gpa: std.mem.Allocator) void {
    var a = pw.Aligner(pw.FnMatch(zbio.scores.blosum62)).init(gpa, -10, -1, .{});
    defer a.deinit();
    var aln = a.local("LSPADKTNVKAA", "PEEKSAV") catch |e| std.debug.panic("blosum: {t}", .{e});
    defer aln.deinit(gpa);
    std.mem.doNotOptimizeAway(aln.score);
}

fn benchHamming(_: std.mem.Allocator) void {
    std.mem.doNotOptimizeAway(zbio.alignment.distance.hamming(&seq_a, &seq_b));
}

fn benchLevenshtein(gpa: std.mem.Allocator) void {
    const d = zbio.alignment.distance.levenshtein(gpa, seq_a[0..512], seq_b[0..512]) catch |e|
        std.debug.panic("levenshtein: {t}", .{e});
    std.mem.doNotOptimizeAway(d);
}

fn benchGcContent(_: std.mem.Allocator) void {
    std.mem.doNotOptimizeAway(zbio.stats.gcContent(&seq_a));
}

fn benchQgrams(_: std.mem.Allocator) void {
    var g = dna_ranks.qgrams(8, &seq_a);
    var acc: usize = 0;
    while (g.next()) |q| acc +%= q;
    std.mem.doNotOptimizeAway(acc);
}

fn fillSeq(buf: []u8, seed: u64) void {
    var prng = std.Random.DefaultPrng.init(seed);
    const rnd = prng.random();
    const bases = "ACGT";
    for (buf) |*c| c.* = bases[rnd.intRangeLessThan(usize, 0, 4)];
}
