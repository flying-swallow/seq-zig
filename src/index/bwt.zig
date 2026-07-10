//! Burrows-Wheeler transform and supporting structures (`less`, `Occ`), ported
//! from rust-bio's `data_structures::bwt`. The alphabet is derived from the BWT
//! itself (its symbol set), so these take just the BWT rather than a separate
//! Alphabet value.

const std = @import("std");
const testing = std.testing;
const simd = @import("../simd.zig");

/// BWT of `text` (which must end in a unique smallest sentinel) given its suffix
/// array `pos`. Owned; caller frees. Complexity O(n).
pub fn bwt(gpa: std.mem.Allocator, text: []const u8, pos: []const usize) ![]u8 {
    std.debug.assert(text.len == pos.len);
    const n = text.len;
    const out = try gpa.alloc(u8, n);
    errdefer gpa.free(out);
    for (0..n) |r| {
        const p = pos[r];
        out[r] = if (p > 0) text[p - 1] else text[n - 1];
    }
    return out;
}

fn maxSymbol(seq: []const u8) u8 {
    var m: u8 = 0;
    for (seq) |c| m = @max(m, c);
    return m;
}

/// The `less` array: `less[c]` = number of BWT symbols lexicographically smaller
/// than `c` (an exclusive prefix-sum of symbol counts). Length maxsym+2.
pub fn less(gpa: std.mem.Allocator, bwt_slice: []const u8) ![]usize {
    const m = @as(usize, maxSymbol(bwt_slice)) + 2;
    const out = try gpa.alloc(usize, m);
    errdefer gpa.free(out);
    @memset(out, 0);
    // Single-pass histogram: already bandwidth-optimal and build-time only. A
    // per-symbol SIMD count would reread the BWT once per symbol, so it is left
    // scalar (SIMD helps the query path via `countByte`, not this build pass).
    for (bwt_slice) |c| out[c] += 1;
    // exclusive prefix scan (rust-bio's `prescan`)
    var s: usize = 0;
    for (out) |*v| {
        const t = v.*;
        v.* = s;
        s += t;
    }
    return out;
}

fn bwtfind(gpa: std.mem.Allocator, bwt_slice: []const u8) ![]usize {
    const n = bwt_slice.len;
    var l = try less(gpa, bwt_slice);
    defer gpa.free(l);
    const find = try gpa.alloc(usize, n);
    errdefer gpa.free(find);
    for (bwt_slice, 0..) |c, r| {
        find[l[c]] = r;
        l[c] += 1;
    }
    return find;
}

/// Invert a BWT back to the original text. Requires the sentinel to be unique
/// and lexicographically smallest.
pub fn invertBwt(gpa: std.mem.Allocator, bwt_slice: []const u8) ![]u8 {
    const n = bwt_slice.len;
    const find = try bwtfind(gpa, bwt_slice);
    defer gpa.free(find);
    const out = try gpa.alloc(u8, n);
    errdefer gpa.free(out);
    var r = find[0];
    for (0..n) |i| {
        r = find[r];
        out[i] = bwt_slice[r];
    }
    return out;
}

fn countByte(slice: []const u8, a: u8) usize {
    // SIMD when the slice spans a full vector; `Occ.get` calls this on slices
    // shorter than the sampling interval `k`, which take the scalar path inside.
    return simd.countEqualByte(slice, a);
}

/// Sampled occurrence array: `get(bwt, r, a)` = count of symbol `a` in
/// `bwt[..r+1]`. Every k-th prefix count is stored. Rows for symbols not present
/// in the BWT are empty.
pub const Occ = struct {
    occ: [][]usize, // rows indexed by symbol byte, length maxsym+1
    k: u32,

    pub fn init(gpa: std.mem.Allocator, bwt_slice: []const u8, k: u32) !Occ {
        const n = bwt_slice.len;
        std.debug.assert(n > 0 and k > 0);
        const m = @as(usize, maxSymbol(bwt_slice)) + 1;

        var present: [256]bool = @splat(false);
        for (bwt_slice) |c| present[c] = true;

        const checkpoints = (n - 1) / k + 1; // count of i in [0,n) with i%k==0

        const occ = try gpa.alloc([]usize, m);
        // Initialize rows so partial cleanup on error is well-defined.
        for (occ) |*row| row.* = &.{};
        errdefer {
            for (occ) |row| if (row.len != 0) gpa.free(row);
            gpa.free(occ);
        }
        for (0..m) |a| {
            if (present[a]) occ[a] = try gpa.alloc(usize, checkpoints);
        }

        const curr = try gpa.alloc(usize, m);
        defer gpa.free(curr);
        @memset(curr, 0);

        // Single-pass running histogram with checkpoints every k-th position.
        // Left scalar: a SIMD reformulation would count each present symbol over
        // tiny (k-sized) blocks, whose per-block reduction overhead outweighs the
        // one-increment-per-byte scalar pass.
        var ci: usize = 0;
        for (bwt_slice, 0..) |c, i| {
            curr[c] += 1;
            if (i % k == 0) {
                for (0..m) |a| {
                    if (present[a]) occ[a][ci] = curr[a];
                }
                ci += 1;
            }
        }

        return .{ .occ = occ, .k = k };
    }

    pub fn deinit(self: *Occ, gpa: std.mem.Allocator) void {
        for (self.occ) |row| if (row.len != 0) gpa.free(row);
        gpa.free(self.occ);
        self.* = undefined;
    }

    pub fn get(self: *const Occ, bwt_slice: []const u8, r: usize, a: u8) usize {
        const k: usize = self.k;
        const lo_checkpoint = r / k;
        const row = self.occ[a];
        const lo_occ = row[lo_checkpoint];

        if (self.k > 64) {
            const hi_checkpoint = lo_checkpoint + 1;
            if (hi_checkpoint < row.len) {
                const hi_occ = row[hi_checkpoint];
                if (lo_occ == hi_occ) return lo_occ;
                const hi_idx = hi_checkpoint * k;
                if ((hi_idx - r) < (k / 2)) {
                    return hi_occ - countByte(bwt_slice[r + 1 .. hi_idx + 1], a);
                }
            }
        }

        const lo_idx = lo_checkpoint * k;
        return countByte(bwt_slice[lo_idx + 1 .. r + 1], a) + lo_occ;
    }
};

// Tests -----------------------------------------------------------------------

const suffix_array = @import("suffix_array.zig");

test "bwt canonical" {
    const text = "GCCTTAACATTATTACGCCTA$";
    const pos = try suffix_array.suffixArray(testing.allocator, text);
    defer testing.allocator.free(pos);
    const b = try bwt(testing.allocator, text, pos);
    defer testing.allocator.free(b);
    try testing.expectEqualStrings("ATTATTCAGGACCC$CTTTCAA", b);
}

test "bwtfind and invert round trip" {
    const text = "cabca$";
    const pos = try suffix_array.suffixArray(testing.allocator, text);
    defer testing.allocator.free(pos);
    const b = try bwt(testing.allocator, text, pos);
    defer testing.allocator.free(b);

    const find = try bwtfind(testing.allocator, b);
    defer testing.allocator.free(find);
    try testing.expectEqualSlices(usize, &.{ 5, 0, 3, 4, 1, 2 }, find);

    const inv = try invertBwt(testing.allocator, b);
    defer testing.allocator.free(inv);
    try testing.expectEqualStrings(text, inv);
}

test "occ sampled table and get" {
    const b = [_]u8{ 1, 3, 3, 1, 2, 0 };
    var occ = try Occ.init(testing.allocator, &b, 3);
    defer occ.deinit(testing.allocator);
    try testing.expectEqualSlices(usize, &.{ 0, 0 }, occ.occ[0]);
    try testing.expectEqualSlices(usize, &.{ 1, 2 }, occ.occ[1]);
    try testing.expectEqualSlices(usize, &.{ 0, 0 }, occ.occ[2]);
    try testing.expectEqualSlices(usize, &.{ 0, 2 }, occ.occ[3]);
    try testing.expectEqual(@as(usize, 1), occ.get(&b, 4, 2));
    try testing.expectEqual(@as(usize, 2), occ.get(&b, 4, 3));
}

test "occ get matches brute force (k>64 branch)" {
    const text = "GCCTTAACATTATTACGCCTA$";
    const pos = try suffix_array.suffixArray(testing.allocator, text);
    defer testing.allocator.free(pos);
    const b = try bwt(testing.allocator, text, pos);
    defer testing.allocator.free(b);
    // Use a large sampling rate to exercise the k>64 backward-count branch.
    var occ = try Occ.init(testing.allocator, b, 128);
    defer occ.deinit(testing.allocator);
    for ("ACGT$") |c| {
        for (0..b.len) |p| {
            const brute = countByte(b[0 .. p + 1], c);
            try testing.expectEqual(brute, occ.get(b, p, c));
        }
    }
}
