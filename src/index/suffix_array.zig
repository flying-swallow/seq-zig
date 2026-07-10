//! Suffix array and LCP array.
//!
//! Deviation from the plan: rust-bio constructs the suffix array with SAIS
//! (linear time). Here it is built by prefix doubling (Manber-Myers, O(n log^2 n)).
//! The suffix array of a text is *unique*, so this yields byte-identical results
//! to rust-bio's SAIS on every test vector, while sidestepping SAIS's
//! wrapping-arithmetic/sentinel hazards (the plan's highest-risk item). SAIS
//! could replace this later for O(n) construction without changing the API.
//!
//! The LCP array uses Kasai's algorithm and is returned uncompressed as []isize.

const std = @import("std");
const testing = std.testing;

/// Suffix array of a byte text. The text should end in a unique, lexicographically
/// smallest sentinel (e.g. '$'). Returns owned positions; caller frees.
pub fn suffixArray(gpa: std.mem.Allocator, text: []const u8) ![]usize {
    return suffixArrayInt(u8, gpa, text);
}

/// Suffix array over a slice of unsigned integers (rust-bio's `suffix_array_int`).
pub fn suffixArrayInt(comptime T: type, gpa: std.mem.Allocator, text: []const T) ![]usize {
    const n = text.len;
    const sa = try gpa.alloc(usize, n);
    errdefer gpa.free(sa);
    if (n == 0) return sa;

    var rank = try gpa.alloc(usize, n);
    defer gpa.free(rank);
    var tmp = try gpa.alloc(usize, n);
    defer gpa.free(tmp);

    for (0..n) |i| {
        sa[i] = i;
        rank[i] = @intCast(text[i]);
    }

    var k: usize = 1;
    while (true) {
        const ctx = Ctx{ .rank = rank, .k = k, .n = n };
        std.mem.sort(usize, sa, ctx, saLess);

        // Recompute ranks from the freshly sorted order.
        tmp[sa[0]] = 0;
        for (1..n) |i| {
            tmp[sa[i]] = tmp[sa[i - 1]] + @intFromBool(saLess(ctx, sa[i - 1], sa[i]));
        }
        @memcpy(rank, tmp);

        if (rank[sa[n - 1]] == n - 1) break; // all suffixes distinct
        k *= 2;
        if (k >= n) break;
    }
    return sa;
}

const Ctx = struct { rank: []const usize, k: usize, n: usize };

fn secondKey(ctx: Ctx, i: usize) isize {
    return if (i + ctx.k < ctx.n) @intCast(ctx.rank[i + ctx.k]) else -1;
}

fn saLess(ctx: Ctx, a: usize, b: usize) bool {
    if (ctx.rank[a] != ctx.rank[b]) return ctx.rank[a] < ctx.rank[b];
    return secondKey(ctx, a) < secondKey(ctx, b);
}

/// LCP array via Kasai's algorithm, uncompressed. Length n+1; index 0 and index
/// n stay -1 (the sentinel suffix and the trailing pad). `lcp[rank]` is the
/// longest common prefix of the suffix at SA rank `rank` and its predecessor.
pub fn lcp(gpa: std.mem.Allocator, text: []const u8, pos: []const usize) ![]isize {
    std.debug.assert(text.len == pos.len);
    const n = text.len;

    var rank = try gpa.alloc(usize, n);
    defer gpa.free(rank);
    for (pos, 0..) |p, r| rank[p] = r;

    const out = try gpa.alloc(isize, n + 1);
    errdefer gpa.free(out);
    @memset(out, -1);

    var h: usize = 0;
    for (0..n) |p| {
        const r = rank[p];
        if (r > 0) {
            const pred = pos[r - 1];
            while (p + h < n and pred + h < n and text[p + h] == text[pred + h]) h += 1;
            out[r] = @intCast(h);
            if (h > 0) h -= 1;
        } else {
            h = 0;
        }
    }
    return out;
}

// Tests -----------------------------------------------------------------------

test "suffix array canonical" {
    const pos = try suffixArray(testing.allocator, "GCCTTAACATTATTACGCCTA$");
    defer testing.allocator.free(pos);
    try testing.expectEqualSlices(usize, &.{
        21, 20, 5, 6, 14, 11, 8, 7, 17, 1, 15, 18, 2, 16, 0, 19, 4, 13, 10, 3, 12, 9,
    }, pos);
}

test "suffix array issue10" {
    {
        const pos = try suffixArray(testing.allocator, "TGTGTGTGTG$");
        defer testing.allocator.free(pos);
        try testing.expectEqualSlices(usize, &.{ 10, 9, 7, 5, 3, 1, 8, 6, 4, 2, 0 }, pos);
    }
    {
        const pos = try suffixArray(testing.allocator, "TGTGTGTG$");
        defer testing.allocator.free(pos);
        try testing.expectEqualSlices(usize, &.{ 8, 7, 5, 3, 1, 6, 4, 2, 0 }, pos);
    }
}

test "suffix array int" {
    const text = [_]usize{ 3, 2, 2, 4, 4, 1, 2, 1, 0 };
    const pos = try suffixArrayInt(usize, testing.allocator, &text);
    defer testing.allocator.free(pos);
    try testing.expectEqualSlices(usize, &.{ 8, 7, 5, 6, 1, 2, 0, 4, 3 }, pos);
}

test "suffix array sorts suffixes lexicographically (property)" {
    // Random DNA-with-sentinel; assert SA order is lexicographic.
    var prng = std.Random.DefaultPrng.init(0xBEEF);
    const rnd = prng.random();
    const bases = "ACGT";
    var trial: usize = 0;
    while (trial < 20) : (trial += 1) {
        const len = rnd.intRangeAtMost(usize, 1, 200);
        const buf = try testing.allocator.alloc(u8, len + 1);
        defer testing.allocator.free(buf);
        for (buf[0..len]) |*c| c.* = bases[rnd.intRangeLessThan(usize, 0, 4)];
        buf[len] = '$';
        const pos = try suffixArray(testing.allocator, buf);
        defer testing.allocator.free(pos);
        for (1..pos.len) |i| {
            try testing.expect(std.mem.lessThan(u8, buf[pos[i - 1]..], buf[pos[i]..]));
        }
    }
}

test "lcp canonical" {
    const text = "GCCTTAACATTATTACGCCTA$";
    const pos = try suffixArray(testing.allocator, text);
    defer testing.allocator.free(pos);
    const l = try lcp(testing.allocator, text, pos);
    defer testing.allocator.free(l);
    try testing.expectEqualSlices(isize, &.{
        -1, 0, 1, 1, 2, 1, 4, 0, 1, 3, 1, 1, 2, 0, 4, 0, 2, 2, 2, 1, 3, 3, -1,
    }, l);
    try testing.expectEqual(@as(isize, 4), l[6]);
}
