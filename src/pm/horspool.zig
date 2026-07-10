//! Horspool exact pattern matching (a simpler, often-faster Boyer-Moore).
//! Ported from rust-bio's `pattern_matching::horspool`. Best case O(n/m), worst
//! O(n*m). No allocation — the shift table is a fixed [256]usize.

const std = @import("std");
const testing = std.testing;

pub const Horspool = struct {
    shift: [256]usize,
    m: usize,
    pattern: []const u8,

    /// Build a matcher for `pattern` (must be non-empty).
    pub fn init(pattern: []const u8) Horspool {
        std.debug.assert(pattern.len > 0);
        const m = pattern.len;
        var shift: [256]usize = @splat(m); // m for symbols not in pattern[..m-1]
        for (pattern[0 .. m - 1], 0..) |a, j| shift[a] = m - 1 - j;
        return .{ .shift = shift, .m = m, .pattern = pattern };
    }

    pub fn findAll(self: *const Horspool, text: []const u8) Matches {
        return .{
            .horspool = self,
            .text = text,
            .n = text.len,
            .last = self.m - 1,
            .pattern_last = self.pattern[self.m - 1],
        };
    }

    pub const Matches = struct {
        horspool: *const Horspool,
        text: []const u8,
        n: usize,
        last: usize,
        pattern_last: u8,

        pub fn next(self: *Matches) ?usize {
            const h = self.horspool;
            while (true) {
                // Shift until the last window symbol matches the pattern's last.
                while (self.last < self.n and self.text[self.last] != self.pattern_last) {
                    self.last += h.shift[self.text[self.last]];
                }
                if (self.last >= self.n) return null;

                const i = self.last + 1 - h.m; // putative start
                const j = self.last;
                self.last += h.shift[self.pattern_last]; // shift again

                if (std.mem.eql(u8, self.text[i..j], h.pattern[0 .. h.m - 1])) return i;
            }
        }
    };
};

// Tests -----------------------------------------------------------------------

fn collect(gpa: std.mem.Allocator, hs: *const Horspool, text: []const u8) ![]usize {
    var out: std.ArrayListUnmanaged(usize) = .empty;
    var it = hs.findAll(text);
    while (it.next()) |p| try out.append(gpa, p);
    return out.toOwnedSlice(gpa);
}

test "horspool shift table" {
    const hs = Horspool.init("AACB");
    try testing.expectEqual(@as(usize, 2), hs.shift['A']);
    try testing.expectEqual(@as(usize, 1), hs.shift['C']);
    try testing.expectEqual(@as(usize, 4), hs.shift['B']);
    try testing.expectEqual(@as(usize, 4), hs.shift['X']);
}

test "horspool find_all" {
    {
        const hs = Horspool.init("qnnnannan");
        const occ = try collect(testing.allocator, &hs, "dhjalkjwqnnnannanaflkjdklfj");
        defer testing.allocator.free(occ);
        try testing.expectEqualSlices(usize, &.{8}, occ);
    }
    {
        const hs = Horspool.init("dhjalk");
        const occ = try collect(testing.allocator, &hs, "dhjalkjwqnnnannanaflkjdklfj");
        defer testing.allocator.free(occ);
        try testing.expectEqualSlices(usize, &.{0}, occ);
    }
    {
        const hs = Horspool.init("GAAAA");
        const occ = try collect(testing.allocator, &hs, "ACGGCTAGGAAAAAGACTGAGGACTGAAAA");
        defer testing.allocator.free(occ);
        try testing.expectEqualSlices(usize, &.{ 8, 25 }, occ);
    }
}
