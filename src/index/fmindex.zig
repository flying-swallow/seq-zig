//! FM-index for finding suffix-array intervals matching a pattern in O(m) time.
//! Ported from rust-bio's `data_structures::fmindex` (the `FMIndex` +
//! `backward_search`). The bidirectional FMD-index (SMEMs over DNA strands) is
//! deferred, as is Myers approximate matching.

const std = @import("std");
const testing = std.testing;
const bwtmod = @import("bwt.zig");

/// A suffix-array interval [lower, upper).
pub const Interval = struct {
    lower: usize,
    upper: usize,

    /// Materialize the matched text positions via the suffix array. Owned.
    pub fn occ(self: Interval, gpa: std.mem.Allocator, sa: []const usize) ![]usize {
        return gpa.dupe(usize, sa[self.lower..self.upper]);
    }
};

/// Result of a backward search.
pub const BackwardSearchResult = union(enum) {
    /// The whole pattern matched; interval is its SA range.
    complete: Interval,
    /// Only a maximal suffix of length `len` matched; interval is its SA range.
    partial: struct { interval: Interval, len: usize },
    /// No suffix of the pattern matched.
    absent,
};

/// FM-index over a BWT plus its `less` and `Occ` structures (borrowed).
pub const FMIndex = struct {
    bwt: []const u8,
    less: []const usize,
    occ: *const bwtmod.Occ,

    pub fn init(bwt_slice: []const u8, less_arr: []const usize, occ: *const bwtmod.Occ) FMIndex {
        return .{ .bwt = bwt_slice, .less = less_arr, .occ = occ };
    }

    fn occAt(self: *const FMIndex, r: usize, a: u8) usize {
        return self.occ.get(self.bwt, r, a);
    }
    fn lessOf(self: *const FMIndex, a: u8) usize {
        return self.less[a];
    }

    /// Backward search for `pattern`. Complexity O(m). See `BackwardSearchResult`.
    pub fn backwardSearch(self: *const FMIndex, pattern: []const u8) BackwardSearchResult {
        var l: usize = 0;
        var r: usize = self.bwt.len - 1;
        var pl = l;
        var pr = r;
        var matched_len: usize = 0;
        var complete_match = true;

        var idx = pattern.len;
        while (idx > 0) {
            idx -= 1;
            const a = pattern[idx];
            const less = self.lessOf(a);
            pl = l;
            pr = r;
            const occ_r = self.occAt(r, a);
            // Empty interval; the `r` assignment would underflow when less==0.
            // (rust-bio issue #606 — the guard is load-bearing in Zig too.)
            if (occ_r == 0) {
                complete_match = false;
                break;
            }
            l = less + (if (l > 0) self.occAt(l - 1, a) else 0);
            r = less + occ_r - 1;
            if (l > r) {
                complete_match = false;
                break;
            }
            matched_len += 1;
        }

        if (matched_len > 0) {
            if (complete_match) {
                return .{ .complete = .{ .lower = l, .upper = r + 1 } };
            }
            return .{ .partial = .{ .interval = .{ .lower = pl, .upper = pr + 1 }, .len = matched_len } };
        }
        return .absent;
    }
};

// Tests -----------------------------------------------------------------------

const suffix_array = @import("suffix_array.zig");

const Fm = struct {
    sa: []usize,
    bwt: []u8,
    less: []usize,
    occ: bwtmod.Occ,
    fm: FMIndex,

    fn build(gpa: std.mem.Allocator, text: []const u8, k: u32) !Fm {
        const sa = try suffix_array.suffixArray(gpa, text);
        const b = try bwtmod.bwt(gpa, text, sa);
        const l = try bwtmod.less(gpa, b);
        const occ = try bwtmod.Occ.init(gpa, b, k);
        return .{ .sa = sa, .bwt = b, .less = l, .occ = occ, .fm = undefined };
    }
    fn index(self: *Fm) void {
        self.fm = FMIndex.init(self.bwt, self.less, &self.occ);
    }
    fn deinit(self: *Fm, gpa: std.mem.Allocator) void {
        gpa.free(self.sa);
        gpa.free(self.bwt);
        gpa.free(self.less);
        self.occ.deinit(gpa);
    }
};

test "fmindex backward_search complete" {
    var x = try Fm.build(testing.allocator, "GCCTTAACATTATTACGCCTA$", 3);
    defer x.deinit(testing.allocator);
    x.index();

    const bsr = x.fm.backwardSearch("TTA");
    const positions = switch (bsr) {
        .complete => |iv| try iv.occ(testing.allocator, x.sa),
        .partial => |p| try p.interval.occ(testing.allocator, x.sa),
        .absent => try testing.allocator.dupe(usize, &.{}),
    };
    defer testing.allocator.free(positions);

    // Order is SA-interval order: [3, 12, 9].
    try testing.expectEqualSlices(usize, &.{ 3, 12, 9 }, positions);
}

test "fmindex backward_search absent and partial" {
    var x = try Fm.build(testing.allocator, "GCCTTAACATTATTACGCCTA$", 3);
    defer x.deinit(testing.allocator);
    x.index();

    // Empty pattern matches nothing -> Absent (the only clean Absent path; an
    // out-of-alphabet symbol would index `less`/`occ` out of bounds, as in rust).
    try testing.expect(x.fm.backwardSearch("") == .absent);

    // "AAA" does not occur (no triple-A), but the suffix "AA" does -> Partial.
    const bsr = x.fm.backwardSearch("AAA");
    try testing.expect(bsr == .partial);
    try testing.expectEqual(@as(usize, 2), bsr.partial.len);
}
