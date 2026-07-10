//! Internal SIMD helpers for the data-parallel byte kernels (Hamming distance,
//! FM-index occurrence counts, GC content). Each kernel is portable: the vector
//! width comes from `std.simd.suggestVectorLength`, which returns `null` on
//! targets without SIMD — there we fall back to a plain scalar loop, so results
//! are identical everywhere. Unlike rust-bio's `distance::simd` (which needs the
//! `triple_accel` crate and runtime CPU detection), `@Vector` lowers to whatever
//! the compile target supports, so this stays dependency-free.
//!
//! Not part of the public API — used by `alignment/distance`, `index/bwt`, and
//! `stats`.

const std = @import("std");
const testing = std.testing;

/// Vector width for byte kernels, or 0 when the target has no SIMD (scalar path).
const L: usize = std.simd.suggestVectorLength(u8) orelse 0;

/// Horizontal sum of a small lane counter. Called only at flush boundaries, so
/// the unrolled scalar add is negligible and avoids vector-widening subtleties.
fn hsum(comptime n: usize, v: @Vector(n, u8)) usize {
    var s: usize = 0;
    inline for (0..n) |i| s += v[i];
    return s;
}

/// Count bytes of `hay` equal to `needle`. Equivalent to
/// `for (hay) |c| n += @intFromBool(c == needle)`.
pub fn countEqualByte(hay: []const u8, needle: u8) usize {
    if (L == 0 or hay.len < L) { // sub-vector inputs: pure scalar, no SIMD overhead
        var n: usize = 0;
        for (hay) |c| n += @intFromBool(c == needle);
        return n;
    }
    const V = @Vector(L, u8);
    const needle_v: V = @splat(needle);
    const ones: V = @splat(1);
    const zeros: V = @splat(0);

    var total: usize = 0;
    var counts: V = zeros;
    var block: usize = 0;
    var i: usize = 0;
    while (i + L <= hay.len) : (i += L) {
        const chunk: V = hay[i..][0..L].*;
        counts += @select(u8, chunk == needle_v, ones, zeros);
        block += 1;
        if (block == 255) { // a u8 lane holds at most 255 before it must flush
            total += hsum(L, counts);
            counts = zeros;
            block = 0;
        }
    }
    total += hsum(L, counts);
    while (i < hay.len) : (i += 1) total += @intFromBool(hay[i] == needle);
    return total;
}

/// Count positions where `a[i] != b[i]` (Hamming distance). Asserts equal length,
/// matching the `for (a, b)` contract of the scalar version.
pub fn countMismatches(a: []const u8, b: []const u8) usize {
    std.debug.assert(a.len == b.len);
    if (L == 0 or a.len < L) {
        var n: usize = 0;
        for (a, b) |x, y| n += @intFromBool(x != y);
        return n;
    }
    const V = @Vector(L, u8);
    const ones: V = @splat(1);
    const zeros: V = @splat(0);

    var total: usize = 0;
    var counts: V = zeros;
    var block: usize = 0;
    var i: usize = 0;
    while (i + L <= a.len) : (i += L) {
        const av: V = a[i..][0..L].*;
        const bv: V = b[i..][0..L].*;
        counts += @select(u8, av != bv, ones, zeros);
        block += 1;
        if (block == 255) {
            total += hsum(L, counts);
            counts = zeros;
            block = 0;
        }
    }
    total += hsum(L, counts);
    while (i < a.len) : (i += 1) total += @intFromBool(a[i] != b[i]);
    return total;
}

/// Count G/C bases (both cases) in `seq`. Equivalent to the scalar
/// `switch (c) { 'g','G','c','C' => n += 1, else => {} }`: fold case with `| 0x20`
/// (upper 'G'|0x20 == 'g', 'C'|0x20 == 'c'; non-letters can't collide), then match
/// 'g' or 'c'. The two matches are mutually exclusive, so two masks never
/// double-count the same byte.
pub fn countGcBytes(seq: []const u8) usize {
    if (L == 0 or seq.len < L) {
        var n: usize = 0;
        for (seq) |c| n += @intFromBool(isGc(c));
        return n;
    }
    const V = @Vector(L, u8);
    const lower: V = @splat(0x20);
    const g_v: V = @splat('g');
    const c_v: V = @splat('c');
    const ones: V = @splat(1);
    const zeros: V = @splat(0);

    var total: usize = 0;
    var counts: V = zeros;
    var block: usize = 0;
    var i: usize = 0;
    while (i + L <= seq.len) : (i += L) {
        const chunk: V = seq[i..][0..L].*;
        const folded = chunk | lower;
        counts += @select(u8, folded == g_v, ones, zeros);
        counts += @select(u8, folded == c_v, ones, zeros);
        block += 1;
        if (block == 127) { // two adds per block, so cap at 127*2 = 254
            total += hsum(L, counts);
            counts = zeros;
            block = 0;
        }
    }
    total += hsum(L, counts);
    while (i < seq.len) : (i += 1) total += @intFromBool(isGc(seq[i]));
    return total;
}

inline fn isGc(c: u8) bool {
    return switch (c) {
        'g', 'G', 'c', 'C' => true,
        else => false,
    };
}

// ---- tests ------------------------------------------------------------------

fn refEqual(hay: []const u8, needle: u8) usize {
    var n: usize = 0;
    for (hay) |c| n += @intFromBool(c == needle);
    return n;
}
fn refMismatch(a: []const u8, b: []const u8) usize {
    var n: usize = 0;
    for (a, b) |x, y| n += @intFromBool(x != y);
    return n;
}
fn refGc(seq: []const u8) usize {
    var n: usize = 0;
    for (seq) |c| n += @intFromBool(isGc(c));
    return n;
}

test "countEqualByte matches scalar across edge sizes" {
    var prng = std.Random.DefaultPrng.init(0xC0FFEE);
    const rng = prng.random();
    // Cover empty, < vlen, == vlen, non-multiple, and a length past the u8-lane
    // flush boundary (255 * up-to-64-wide vectors).
    const sizes = [_]usize{ 0, 1, 7, 16, 31, 32, 33, 64, 1000, 255 * 64 + 5 };
    var buf: [255 * 64 + 5]u8 = undefined;
    for (sizes) |n| {
        const seq = buf[0..n];
        for (seq) |*c| c.* = "ACGTN"[rng.intRangeLessThan(usize, 0, 5)];
        for ("ACGTNX") |needle| {
            try testing.expectEqual(refEqual(seq, needle), countEqualByte(seq, needle));
        }
    }
}

test "countMismatches matches scalar and handles all/none equal" {
    var prng = std.Random.DefaultPrng.init(0xBEEF);
    const rng = prng.random();
    const sizes = [_]usize{ 0, 1, 15, 32, 100, 255 * 64 + 3 };
    var a: [255 * 64 + 3]u8 = undefined;
    var b: [255 * 64 + 3]u8 = undefined;
    for (sizes) |n| {
        const av = a[0..n];
        const bv = b[0..n];
        for (av, bv) |*x, *y| {
            x.* = "ACGT"[rng.intRangeLessThan(usize, 0, 4)];
            y.* = "ACGT"[rng.intRangeLessThan(usize, 0, 4)];
        }
        try testing.expectEqual(refMismatch(av, bv), countMismatches(av, bv));
        // all equal → 0
        @memcpy(bv, av);
        try testing.expectEqual(@as(usize, 0), countMismatches(av, bv));
        // all different (complement-ish) → n
        for (av, bv) |x, *y| y.* = if (x == 'A') 'C' else 'A';
        try testing.expectEqual(n, countMismatches(av, bv));
    }
}

test "countGcBytes matches scalar including mixed case and flush boundary" {
    var prng = std.Random.DefaultPrng.init(0x5EED);
    const rng = prng.random();
    const sizes = [_]usize{ 0, 1, 8, 32, 63, 500, 127 * 64 + 9 };
    var buf: [127 * 64 + 9]u8 = undefined;
    for (sizes) |n| {
        const seq = buf[0..n];
        for (seq) |*c| c.* = "aAcCgGtTnN-"[rng.intRangeLessThan(usize, 0, 11)];
        try testing.expectEqual(refGc(seq), countGcBytes(seq));
    }
    // all G/C → full length
    var full: [200]u8 = undefined;
    for (&full, 0..) |*c, i| c.* = "gGcC"[i % 4];
    try testing.expectEqual(@as(usize, 200), countGcBytes(&full));
}
