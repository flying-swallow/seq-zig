//! Knuth-Morris-Pratt exact matching. Ported from rust-bio's
//! `pattern_matching::kmp`. Builds the `lps` (longest-proper-prefix-suffix)
//! table and scans the text in O(n). Owns the `lps` heap array — call `deinit`.

const std = @import("std");
const testing = std.testing;

pub const KMP = struct {
    m: usize,
    lps: []usize,
    pattern: []const u8,

    pub fn init(gpa: std.mem.Allocator, pattern: []const u8) !KMP {
        std.debug.assert(pattern.len > 0);
        const lps_tab = try computeLps(gpa, pattern);
        return .{ .m = pattern.len, .lps = lps_tab, .pattern = pattern };
    }

    pub fn deinit(self: *KMP, gpa: std.mem.Allocator) void {
        gpa.free(self.lps);
        self.* = undefined;
    }

    pub fn delta(self: *const KMP, q_in: usize, a: u8) usize {
        var q = q_in;
        while (q == self.m or (self.pattern[q] != a and q > 0)) {
            q = self.lps[q - 1];
        }
        if (self.pattern[q] == a) q += 1;
        return q;
    }

    pub fn findAll(self: *const KMP, text: []const u8) Matches {
        return .{ .kmp = self, .q = 0, .text = text, .pos = 0 };
    }

    pub const Matches = struct {
        kmp: *const KMP,
        q: usize,
        text: []const u8,
        pos: usize,

        pub fn next(self: *Matches) ?usize {
            const k = self.kmp;
            while (self.pos < self.text.len) {
                const i = self.pos;
                const c = self.text[i];
                self.pos += 1;
                self.q = k.delta(self.q, c);
                if (self.q == k.m) return 1 + i - k.m;
            }
            return null;
        }
    };
};

fn computeLps(gpa: std.mem.Allocator, pattern: []const u8) ![]usize {
    const m = pattern.len;
    var lps = try gpa.alloc(usize, m);
    errdefer gpa.free(lps);
    @memset(lps, 0);
    var q: usize = 0;
    var i: usize = 1;
    while (i < m) : (i += 1) {
        while (q > 0 and pattern[q] != pattern[i]) q = lps[q - 1];
        if (pattern[q] == pattern[i]) q += 1;
        lps[i] = q;
    }
    return lps;
}

// Tests -----------------------------------------------------------------------

fn collect(gpa: std.mem.Allocator, k: *const KMP, text: []const u8) ![]usize {
    var out: std.ArrayListUnmanaged(usize) = .empty;
    var it = k.findAll(text);
    while (it.next()) |p| try out.append(gpa, p);
    return out.toOwnedSlice(gpa);
}

test "kmp lps" {
    const lps = try computeLps(testing.allocator, "ababaca");
    defer testing.allocator.free(lps);
    try testing.expectEqualSlices(usize, &.{ 0, 0, 1, 2, 3, 0, 1 }, lps);
}

test "kmp delta" {
    var k = try KMP.init(testing.allocator, "abbab");
    defer k.deinit(testing.allocator);
    try testing.expectEqual(@as(usize, 1), k.delta(0, 'a'));
    try testing.expectEqual(@as(usize, 0), k.delta(0, 'b'));
    try testing.expectEqual(@as(usize, 1), k.delta(1, 'a'));
    try testing.expectEqual(@as(usize, 2), k.delta(1, 'b'));
    try testing.expectEqual(@as(usize, 1), k.delta(2, 'a'));
    try testing.expectEqual(@as(usize, 3), k.delta(2, 'b'));
    try testing.expectEqual(@as(usize, 4), k.delta(3, 'a'));
    try testing.expectEqual(@as(usize, 0), k.delta(3, 'b'));
    try testing.expectEqual(@as(usize, 1), k.delta(4, 'a'));
    try testing.expectEqual(@as(usize, 5), k.delta(4, 'b'));
    try testing.expectEqual(@as(usize, 1), k.delta(5, 'a'));
    try testing.expectEqual(@as(usize, 3), k.delta(5, 'b'));
}

test "kmp find_all" {
    {
        var k = try KMP.init(testing.allocator, "abbab");
        defer k.deinit(testing.allocator);
        const occ = try collect(testing.allocator, &k, "aaaaabbabbbbbbbabbab");
        defer testing.allocator.free(occ);
        try testing.expectEqualSlices(usize, &.{ 4, 15 }, occ);
    }
    {
        var k = try KMP.init(testing.allocator, "qnnnannan");
        defer k.deinit(testing.allocator);
        const occ = try collect(testing.allocator, &k, "dhjalkjwqnnnannanaflkjdklfj");
        defer testing.allocator.free(occ);
        try testing.expectEqualSlices(usize, &.{8}, occ);
    }
    {
        var k = try KMP.init(testing.allocator, "dhjalk");
        defer k.deinit(testing.allocator);
        const occ = try collect(testing.allocator, &k, "dhjalkjwqnnnannanaflkjdklfj");
        defer testing.allocator.free(occ);
        try testing.expectEqualSlices(usize, &.{0}, occ);
    }
}
