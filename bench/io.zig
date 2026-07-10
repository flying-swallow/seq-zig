//! Benchmarks for zbio's existing format-IO code: FASTA/FASTQ parsing, .fai
//! scanning, BAM reading, and length statistics. `bytes_per_run`/`items_per_run`
//! make zBench report throughput (MB/s, records/s, elements/s).

const std = @import("std");
const zbench = @import("zbench");
const zbio = @import("zbio");
const fixtures = @import("fixtures");

pub fn register(bench: *zbench.Benchmark) !void {
    try bench.add("fasta parse ce.fa", benchFasta, .{ .bytes_per_run = fixtures.ce_fa.len });
    try bench.add("fastq parse t3.fq", benchFastq, .{ .bytes_per_run = fixtures.t3_fq.len });
    try bench.add("fai index ce.fa", benchFaiIndex, .{ .bytes_per_run = fixtures.ce_fa.len });
    try bench.add("bam read t1.bam", benchBamRead, .{ .bytes_per_run = fixtures.t1_bam.len });
    try bench.add("length stats 100k", benchLengthStats, .{ .items_per_run = 100_000 });
}

// One record per call. FASTA returns `false` at EOF; FASTQ instead raises
// EndOfStream on the trailing call (the leading takeByte isn't guarded) — this
// loop treats both as done. Discarding sinks keep it measuring the parser, not
// output allocation.
fn parseAll(bytes: []const u8) u64 {
    var reader: std.Io.Reader = .fixed(bytes);
    var nb: [512]u8 = undefined;
    var sb: [4096]u8 = undefined;
    var qb: [4096]u8 = undefined;
    var name: std.Io.Writer.Discarding = .init(&nb);
    var seq: std.Io.Writer.Discarding = .init(&sb);
    var qual: std.Io.Writer.Discarding = .init(&qb);

    while (true) {
        const more = zbio.fa.takeFqSequence(&reader, &name.writer, &seq.writer, &qual.writer) catch |e|
            switch (e) {
                error.EndOfStream => break,
                else => std.debug.panic("parse failed: {t}", .{e}),
            };
        if (!more) break;
    }
    return seq.fullCount();
}

fn benchFasta(_: std.mem.Allocator) void {
    std.mem.doNotOptimizeAway(parseAll(fixtures.ce_fa));
}

fn benchFastq(_: std.mem.Allocator) void {
    std.mem.doNotOptimizeAway(parseAll(fixtures.t3_fq));
}

// Scans the FASTA stream and computes the .fai offset/line-width index (the same
// work `samtools faidx` does), not parsing an existing .fai.
fn benchFaiIndex(gpa: std.mem.Allocator) void {
    var reader: std.Io.Reader = .fixed(fixtures.ce_fa);
    var sc = zbio.fai.FaiFaIndexIndexScanner.init(gpa, &reader);
    defer sc.deinit();
    var n: usize = 0;
    while (sc.next() catch |e| std.debug.panic("fai: {t}", .{e})) |_| n +%= 1;
    std.mem.doNotOptimizeAway(n);
}

fn benchBamRead(gpa: std.mem.Allocator) void {
    var reader: std.Io.Reader = .fixed(fixtures.t1_bam);
    var r = zbio.bam.open(gpa, &reader) catch |e| std.debug.panic("bam open: {t}", .{e});
    defer r.deinit();
    var n: usize = 0;
    while (r.next() catch |e| std.debug.panic("bam next: {t}", .{e})) |_| n +%= 1;
    std.mem.doNotOptimizeAway(n);
}

fn benchLengthStats(gpa: std.mem.Allocator) void {
    var s = zbio.stats.LengthStats.init(gpa);
    defer s.deinit();
    var prng = std.Random.DefaultPrng.init(0xC0FFEE);
    const rnd = prng.random();
    var i: usize = 0;
    while (i < 100_000) : (i += 1)
        s.add(rnd.intRangeAtMost(u64, 1, 50_000)) catch |e| std.debug.panic("stats: {t}", .{e});
    std.mem.doNotOptimizeAway(s.n50());
}
