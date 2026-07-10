//! Sequence distance functions, ported from rust-bio's `alignment::distance`.
//! Hamming uses the SIMD byte kernel in `simd.zig` (portable `@Vector`, no
//! dependency). rust-bio delegates Levenshtein to the `editdistancek`/
//! `triple_accel` SIMD crates; to keep zbio dependency-free those stay plain
//! scalar DP.

const std = @import("std");
const testing = std.testing;
const simd = @import("../simd.zig");

/// Hamming distance between two equal-length texts. Complexity O(n).
/// Asserts equal length (matching rust-bio, which panics otherwise).
pub fn hamming(alpha: []const u8, beta: []const u8) u64 {
    std.debug.assert(alpha.len == beta.len);
    return simd.countMismatches(alpha, beta);
}

/// Levenshtein (edit) distance between two texts. Complexity O(n*m) time,
/// O(min(n,m)) space (two rolling rows over the shorter text).
pub fn levenshtein(gpa: std.mem.Allocator, alpha: []const u8, beta: []const u8) !u32 {
    // Iterate columns over the shorter text to bound the row buffers.
    const a, const b = if (alpha.len <= beta.len) .{ alpha, beta } else .{ beta, alpha };
    const m = a.len;

    var prev = try gpa.alloc(u32, m + 1);
    defer gpa.free(prev);
    var curr = try gpa.alloc(u32, m + 1);
    defer gpa.free(curr);

    for (0..m + 1) |j| prev[j] = @intCast(j);

    for (1..b.len + 1) |i| {
        curr[0] = @intCast(i);
        for (1..m + 1) |j| {
            const cost: u32 = if (b[i - 1] == a[j - 1]) 0 else 1;
            curr[j] = @min(@min(prev[j] + 1, curr[j - 1] + 1), prev[j - 1] + cost);
        }
        std.mem.swap([]u32, &prev, &curr);
    }
    return prev[m];
}

/// Bounded Levenshtein distance: returns the edit distance if it is `<= k`,
/// otherwise `null`. Uses the Ukkonen early-out — the minimum value in each DP
/// row is a lower bound on the final distance, so once a whole row exceeds `k`
/// the answer must too.
pub fn boundedLevenshtein(gpa: std.mem.Allocator, alpha: []const u8, beta: []const u8, k: u32) !?u32 {
    const a, const b = if (alpha.len <= beta.len) .{ alpha, beta } else .{ beta, alpha };
    const m = a.len;

    // A length difference greater than k can never be bridged within k edits.
    const len_diff = b.len - m; // b is the longer text
    if (len_diff > k) return null;

    var prev = try gpa.alloc(u32, m + 1);
    defer gpa.free(prev);
    var curr = try gpa.alloc(u32, m + 1);
    defer gpa.free(curr);

    for (0..m + 1) |j| prev[j] = @intCast(j);

    for (1..b.len + 1) |i| {
        curr[0] = @intCast(i);
        var row_min: u32 = curr[0];
        for (1..m + 1) |j| {
            const cost: u32 = if (b[i - 1] == a[j - 1]) 0 else 1;
            curr[j] = @min(@min(prev[j] + 1, curr[j - 1] + 1), prev[j - 1] + cost);
            row_min = @min(row_min, curr[j]);
        }
        std.mem.swap([]u32, &prev, &curr);
        if (row_min > k) return null;
    }
    const dist = prev[m];
    return if (dist > k) null else dist;
}

// Tests -----------------------------------------------------------------------

test "hamming distance" {
    try testing.expectEqual(@as(u64, 5), hamming("GTCTGCATGCG", "TTTAGCTAGCG"));
    try testing.expectEqual(@as(u64, 0), hamming("ACGT", "ACGT"));
}

test "levenshtein distance" {
    try testing.expectEqual(@as(u32, 5), try levenshtein(testing.allocator, "ACCGTGGAT", "AAAAACCGTTGAT"));
    try testing.expectEqual(
        try levenshtein(testing.allocator, "ACCGTGGAT", "AAAAACCGTTGAT"),
        try levenshtein(testing.allocator, "AAAAACCGTTGAT", "ACCGTGGAT"),
    );
    try testing.expectEqual(@as(u32, 4), try levenshtein(testing.allocator, "AAA", "TTTT"));
    try testing.expectEqual(@as(u32, 4), try levenshtein(testing.allocator, "TTTT", "AAA"));
    try testing.expectEqual(@as(u32, 0), try levenshtein(testing.allocator, "", ""));
    try testing.expectEqual(@as(u32, 3), try levenshtein(testing.allocator, "", "ABC"));
}

test "bounded levenshtein distance" {
    const x = "ACCGTGGAT";
    const y = "AAAAACCGTTGAT";
    try testing.expectEqual(@as(?u32, 5), try boundedLevenshtein(testing.allocator, x, y, std.math.maxInt(u32)));
    try testing.expectEqual(@as(?u32, 5), try boundedLevenshtein(testing.allocator, x, y, 5));
    try testing.expectEqual(@as(?u32, null), try boundedLevenshtein(testing.allocator, x, y, 4));
    try testing.expectEqual(@as(?u32, 4), try boundedLevenshtein(testing.allocator, "AAA", "TTTT", 10));
    try testing.expectEqual(@as(?u32, null), try boundedLevenshtein(testing.allocator, "AAA", "TTTTTTTT", 2));
}
