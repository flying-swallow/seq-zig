const std = @import("std");
const testing = std.testing;
const simd = @import("simd.zig");

// Sequence-length statistics and base-composition counting, modeled on
// seqkit's `util.LengthStats` and `byteutil.CountBytes`. These are the
// reusable primitives behind an `seqkit stats`-style command: accumulate a
// length per record, then query count/sum/min/max/mean, quartiles, and the
// assembly metrics N50/L50/NX.

/// Accumulates per-record sequence lengths and computes summary statistics
/// (count/sum/min/max/mean, quartiles, N50/L50, NX).
///
/// Ported to match seqkit's `bio/util.LengthStats` semantics exactly: lengths
/// are grouped into distinct length classes (`[length, count]` bins), and the
/// quartile / N50 metrics are defined over those bins. Notably this makes L50
/// the number of distinct length classes (from longest) needed to reach the N50
/// threshold — not the number of sequences — matching seqkit's output.
pub const LengthStats = struct {
    gpa: std.mem.Allocator,
    lens: std.AutoHashMap(u64, u64), // length -> occurrence count
    total: u64,
    n: u64,
    min_len: u64,
    max_len: u64,

    // Built lazily by `finalize`: distinct lengths ascending, and the running
    // cumulative count alongside each.
    counts: ?[][2]u64,
    acc_counts: ?[][2]u64,
    finalized: bool,

    // N50/L50 memoization.
    l50_val: usize,
    n50_done: bool,

    pub fn init(gpa: std.mem.Allocator) LengthStats {
        return .{
            .gpa = gpa,
            .lens = std.AutoHashMap(u64, u64).init(gpa),
            .total = 0,
            .n = 0,
            .min_len = std.math.maxInt(u64),
            .max_len = 0,
            .counts = null,
            .acc_counts = null,
            .finalized = false,
            .l50_val = 0,
            .n50_done = false,
        };
    }

    pub fn deinit(self: *LengthStats) void {
        self.lens.deinit();
        if (self.counts) |c| self.gpa.free(c);
        if (self.acc_counts) |a| self.gpa.free(a);
    }

    pub fn add(self: *LengthStats, len: u64) !void {
        self.n += 1;
        self.total += len;
        const gop = try self.lens.getOrPut(len);
        if (gop.found_existing) {
            gop.value_ptr.* += 1;
        } else {
            gop.value_ptr.* = 1;
        }
        if (len > self.max_len) self.max_len = len;
        if (len < self.min_len) self.min_len = len;
        self.finalized = false;
        self.n50_done = false;
    }

    pub fn count(self: *const LengthStats) usize {
        return @intCast(self.n);
    }

    pub fn sum(self: *const LengthStats) u64 {
        return self.total;
    }

    pub fn min(self: *const LengthStats) u64 {
        if (self.n == 0) return 0;
        return self.min_len;
    }

    pub fn max(self: *const LengthStats) u64 {
        return self.max_len;
    }

    pub fn mean(self: *const LengthStats) f64 {
        if (self.n == 0) return 0;
        return @as(f64, @floatFromInt(self.total)) / @as(f64, @floatFromInt(self.n));
    }

    fn lessThanBin(_: void, a: [2]u64, b: [2]u64) bool {
        return a[0] < b[0];
    }

    /// Build the sorted distinct-length bins and their cumulative counts.
    fn finalize(self: *LengthStats) !void {
        if (self.finalized) return;
        if (self.counts) |c| self.gpa.free(c);
        if (self.acc_counts) |a| self.gpa.free(a);
        self.counts = null;
        self.acc_counts = null;

        const ndistinct = self.lens.count();
        if (ndistinct == 0) {
            self.finalized = true;
            return;
        }

        const counts = try self.gpa.alloc([2]u64, ndistinct);
        errdefer self.gpa.free(counts);
        var it = self.lens.iterator();
        var i: usize = 0;
        while (it.next()) |entry| : (i += 1) {
            counts[i] = .{ entry.key_ptr.*, entry.value_ptr.* };
        }
        std.sort.pdq([2]u64, counts, {}, lessThanBin);

        const acc = try self.gpa.alloc([2]u64, ndistinct);
        var running: u64 = 0;
        for (counts, 0..) |bin, k| {
            running += bin[1];
            acc[k] = .{ bin[0], running };
        }

        self.counts = counts;
        self.acc_counts = acc;
        self.finalized = true;
    }

    fn getValue(self: *LengthStats, even: bool, i_median_l: u64, i_median_r: u64) f64 {
        const acc = self.acc_counts orelse return 0;
        var flag = false;
        var prev: u64 = 0;
        for (acc) |data| {
            const acc_count = data[1];
            if (flag) {
                // The middle two elements have different lengths; average them.
                return @as(f64, @floatFromInt(data[0] + prev)) / 2.0;
            }
            if (acc_count >= i_median_l + 1) {
                if (even) {
                    if (acc_count >= i_median_r + 1) return @floatFromInt(data[0]);
                    flag = true;
                    prev = data[0];
                } else {
                    return @floatFromInt(data[0]);
                }
            }
        }
        return 0;
    }

    /// Second quartile (median).
    pub fn q2(self: *LengthStats) f64 {
        self.finalize() catch return 0;
        const counts = self.counts orelse return 0;
        if (counts.len == 0) return 0;
        if (counts.len == 1) return @floatFromInt(counts[0][0]);

        const even = self.n & 1 == 0;
        var l: u64 = 0;
        var r: u64 = 0;
        if (even) {
            l = self.n / 2 - 1;
            r = self.n / 2;
        } else {
            l = self.n / 2;
        }
        return self.getValue(even, l, r);
    }

    /// First quartile.
    pub fn q1(self: *LengthStats) f64 {
        self.finalize() catch return 0;
        const counts = self.counts orelse return 0;
        if (counts.len == 0) return 0;
        if (counts.len == 1) return @floatFromInt(counts[0][0]);

        var even = self.n & 1 == 0;
        const nn: u64 = if (even) self.n / 2 else (self.n + 1) / 2;
        even = nn % 2 == 0;
        var l: u64 = 0;
        var r: u64 = 0;
        if (even) {
            l = nn / 2 - 1;
            r = nn / 2;
        } else {
            l = nn / 2;
        }
        return self.getValue(even, l, r);
    }

    /// Third quartile.
    pub fn q3(self: *LengthStats) f64 {
        self.finalize() catch return 0;
        const counts = self.counts orelse return 0;
        if (counts.len == 0) return 0;
        if (counts.len == 1) return @floatFromInt(counts[0][0]);

        var even = self.n & 1 == 0;
        var mean_off: u64 = 0;
        var nn: u64 = 0;
        if (even) {
            nn = self.n / 2;
            mean_off = nn;
        } else {
            nn = (self.n + 1) / 2;
            mean_off = self.n / 2;
        }
        even = nn % 2 == 0;
        var l: u64 = 0;
        var r: u64 = 0;
        if (even) {
            l = nn / 2 - 1 + mean_off;
            r = nn / 2 + mean_off;
        } else {
            l = nn / 2 + mean_off;
        }
        return self.getValue(even, l, r);
    }

    /// N50: the length of the shortest distinct length class among the longest
    /// classes that together cover at least half of the total length. Also
    /// records L50 (number of such distinct classes).
    pub fn n50(self: *LengthStats) u64 {
        self.finalize() catch return 0;
        const counts = self.counts orelse return 0;
        if (counts.len == 0) return 0;
        if (counts.len == 1) {
            self.l50_val = 1;
            self.n50_done = true;
            return counts[0][0];
        }

        const half: f64 = @as(f64, @floatFromInt(self.total)) / 2.0;
        var sum_len: f64 = 0;
        var i: usize = counts.len;
        while (i > 0) {
            i -= 1;
            const data = counts[i];
            sum_len += @floatFromInt(data[0] * data[1]);
            if (sum_len >= half) {
                self.l50_val = counts.len - i;
                self.n50_done = true;
                return data[0];
            }
        }
        return 0;
    }

    /// L50: number of distinct length classes (from longest) needed to reach N50.
    pub fn l50(self: *LengthStats) usize {
        if (!self.n50_done) _ = self.n50();
        return self.l50_val;
    }

    /// NX statistic (N50 generalized): length class covering at least `x`% of
    /// the total when accumulating from the longest. `x` is in [0, 100].
    pub fn nx(self: *LengthStats, x: f64) u64 {
        self.finalize() catch return 0;
        const counts = self.counts orelse return 0;
        if (counts.len == 0) return 0;
        if (counts.len == 1) return counts[0][0];

        const boundary: f64 = @as(f64, @floatFromInt(self.total)) * x / 100.0;
        var sum_len: f64 = 0;
        var i: usize = counts.len;
        while (i > 0) {
            i -= 1;
            const data = counts[i];
            sum_len += @floatFromInt(data[0] * data[1]);
            if (sum_len >= boundary) return data[0];
        }
        return 0;
    }
};

/// Count the bytes of `seq` that appear in `set` (the `byteutil.CountBytes`
/// equivalent used for GC, gap, and ambiguous-base counting).
pub fn countAny(seq: []const u8, set: []const u8) usize {
    var lut: [256]bool = @splat(false);
    for (set) |c| lut[c] = true;
    var n: usize = 0;
    for (seq) |c| {
        if (lut[c]) n += 1;
    }
    return n;
}

/// Count G/C bases (both cases) in `seq`. SIMD-accelerated (see `simd.zig`).
pub fn gcCount(seq: []const u8) usize {
    return simd.countGcBytes(seq);
}

/// GC ratio counter over every `step`-th base (rust-bio `seq_analysis::gc`).
/// The contiguous `step == 1` case uses the SIMD G/C counter; strided steps stay
/// scalar (an every-`step`-th gather is not a contiguous reduction).
fn gcnContent(seq: []const u8, step: usize) f32 {
    if (step == 1) {
        const count = simd.countGcBytes(seq);
        // Empty input yields NaN (0/0), matching rust-bio.
        return @as(f32, @floatFromInt(count)) / @as(f32, @floatFromInt(seq.len));
    }
    var l: usize = 0;
    var count: usize = 0;
    var i: usize = 0;
    while (i < seq.len) : (i += step) {
        l += 1;
        switch (seq[i]) {
            'g', 'G', 'c', 'C' => count += 1,
            else => {},
        }
    }
    // Empty input yields NaN (0/0), matching rust-bio.
    return @as(f32, @floatFromInt(count)) / @as(f32, @floatFromInt(l));
}

/// Ratio of bases that are G or C (both cases). Complexity O(n).
pub fn gcContent(seq: []const u8) f32 {
    return gcnContent(seq, 1);
}

/// Ratio of every 3rd base (positions 0, 3, 6, ...) that is G or C.
pub fn gc3Content(seq: []const u8) f32 {
    return gcnContent(seq, 3);
}

// Tests ---------------------------------------------------------------------

test "LengthStats basic aggregates" {
    var s = LengthStats.init(testing.allocator);
    defer s.deinit();
    try s.add(10);
    try s.add(30);
    try s.add(20);

    try testing.expectEqual(@as(usize, 3), s.count());
    try testing.expectEqual(@as(u64, 60), s.sum());
    try testing.expectEqual(@as(u64, 10), s.min());
    try testing.expectEqual(@as(u64, 30), s.max());
    try testing.expectApproxEqAbs(@as(f64, 20), s.mean(), 1e-9);
}

test "LengthStats empty" {
    var s = LengthStats.init(testing.allocator);
    defer s.deinit();
    try testing.expectEqual(@as(usize, 0), s.count());
    try testing.expectEqual(@as(u64, 0), s.sum());
    try testing.expectEqual(@as(u64, 0), s.n50());
    try testing.expectApproxEqAbs(@as(f64, 0), s.mean(), 1e-9);
}

test "LengthStats N50 / L50 over distinct length classes" {
    // lengths {2,2,2,3,3,4,8,8}, total 32, half 16. Distinct classes ascending:
    // [2x3, 3x2, 4x1, 8x2]. From the top, 8x2 = 16 >= 16 -> N50 = 8, and L50 is
    // the number of distinct classes used (1), matching seqkit's definition.
    var s = LengthStats.init(testing.allocator);
    defer s.deinit();
    for ([_]u64{ 2, 3, 4, 8, 2, 8, 3, 2 }) |x| try s.add(x);
    try testing.expectEqual(@as(u64, 8), s.n50());
    try testing.expectEqual(@as(usize, 1), s.l50());
    try testing.expectEqual(@as(u64, 8), s.nx(50));
}

test "LengthStats single length class" {
    // 8 sequences all length 227: one distinct class -> Q1=Q2=Q3=227, L50=1.
    var s = LengthStats.init(testing.allocator);
    defer s.deinit();
    for (0..8) |_| try s.add(227);
    try testing.expectApproxEqAbs(@as(f64, 227), s.q1(), 1e-9);
    try testing.expectApproxEqAbs(@as(f64, 227), s.q2(), 1e-9);
    try testing.expectApproxEqAbs(@as(f64, 227), s.q3(), 1e-9);
    try testing.expectEqual(@as(u64, 227), s.n50());
    try testing.expectEqual(@as(usize, 1), s.l50());
}

test "LengthStats quartiles match seqkit getValue semantics" {
    var s = LengthStats.init(testing.allocator);
    defer s.deinit();
    // odd count 1..7 -> Q2 = 4; seqkit's exclusive halves give Q1 = 2.5, Q3 = 5.5
    for ([_]u64{ 7, 1, 5, 3, 6, 2, 4 }) |x| try s.add(x);
    try testing.expectApproxEqAbs(@as(f64, 2.5), s.q1(), 1e-9);
    try testing.expectApproxEqAbs(@as(f64, 4), s.q2(), 1e-9);
    try testing.expectApproxEqAbs(@as(f64, 5.5), s.q3(), 1e-9);

    // even count: {10,20,30,40} -> Q1 = 15, Q2 = 25, Q3 = 35
    var e = LengthStats.init(testing.allocator);
    defer e.deinit();
    for ([_]u64{ 40, 10, 30, 20 }) |x| try e.add(x);
    try testing.expectApproxEqAbs(@as(f64, 15), e.q1(), 1e-9);
    try testing.expectApproxEqAbs(@as(f64, 25), e.q2(), 1e-9);
    try testing.expectApproxEqAbs(@as(f64, 35), e.q3(), 1e-9);
}

test "countAny and gcCount" {
    try testing.expectEqual(@as(usize, 4), gcCount("GGCCatAT"));
    try testing.expectEqual(@as(usize, 2), countAny("ACGT-N.n", "Nn"));
    try testing.expectEqual(@as(usize, 2), countAny("A C-GT", "- ")); // one space, one dash
    try testing.expectEqual(@as(usize, 0), countAny("ACGT", "Xx"));
}

test "gcContent / gc3Content" {
    try testing.expectApproxEqAbs(@as(f32, 0.0), gcContent("ATAT"), 1e-6);
    try testing.expectApproxEqAbs(@as(f32, 0.5), gcContent("ATGC"), 1e-6);
    try testing.expectApproxEqAbs(@as(f32, 1.0), gcContent("GCGC"), 1e-6);
    try testing.expectApproxEqAbs(@as(f32, 2.0 / 8.0), gcContent("GATATACA"), 1e-6);
    // step 3 over "GATATACA" -> {G, A, C} -> 2/3
    try testing.expectApproxEqAbs(@as(f32, 2.0 / 3.0), gc3Content("GATATACA"), 1e-6);
}
