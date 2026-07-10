//! Backward Nondeterministic DAWG Matching (BNDM). Ported from rust-bio's
//! `pattern_matching::bndm`. Reuses the Shift-And masks over the reversed
//! pattern. Best case O(n/m), worst O(n*m). Patterns < 64 symbols. No allocation.

const std = @import("std");
const testing = std.testing;
const shift_and = @import("shift_and.zig");

pub const BNDM = struct {
    m: usize,
    table: [256]u64,
    accept: u64,

    pub fn init(pattern: []const u8) BNDM {
        std.debug.assert(pattern.len < 64); // pattern of less than 64 symbols
        // Nondeterministic suffix automaton = Shift-And over the reversed pattern.
        const mk = shift_and.buildMasks(pattern, true);
        return .{ .m = pattern.len, .table = mk.table, .accept = mk.accept };
    }

    pub fn findAll(self: *const BNDM, text: []const u8) Matches {
        return .{ .bndm = self, .window = self.m, .text = text };
    }

    pub const Matches = struct {
        bndm: *const BNDM,
        window: usize,
        text: []const u8,

        pub fn next(self: *Matches) ?usize {
            const b = self.bndm;
            const m_shift: u6 = @intCast(b.m);
            while (self.window <= self.text.len) {
                var occ: ?usize = null;
                var active: u64 = (@as(u64, 1) << m_shift) - 1; // all states active
                var j: usize = 1;
                var lastsuffix: usize = 0;
                while (active != 0) {
                    active &= b.table[self.text[self.window - j]];
                    if (active & b.accept != 0) {
                        if (j == b.m) {
                            occ = self.window - b.m;
                            break;
                        } else {
                            lastsuffix = j; // a prefix of length j matches
                        }
                    }
                    j += 1;
                    active <<= 1;
                }
                self.window += b.m - lastsuffix;
                if (occ) |o| return o;
            }
            return null;
        }
    };
};

// Tests -----------------------------------------------------------------------

fn collect(gpa: std.mem.Allocator, b: *const BNDM, text: []const u8) ![]usize {
    var out: std.ArrayListUnmanaged(usize) = .empty;
    var it = b.findAll(text);
    while (it.next()) |p| try out.append(gpa, p);
    return out.toOwnedSlice(gpa);
}

test "bndm find_all" {
    {
        const b = BNDM.init("GAAAA");
        const occ = try collect(testing.allocator, &b, "ACGGCTAGAAAAGGCTAGAAAA");
        defer testing.allocator.free(occ);
        try testing.expectEqualSlices(usize, &.{ 7, 17 }, occ);
    }
    {
        const b = BNDM.init("qnnnannan");
        const occ = try collect(testing.allocator, &b, "dhjalkjwqnnnannanaflkjdklfj");
        defer testing.allocator.free(occ);
        try testing.expectEqualSlices(usize, &.{8}, occ);
    }
    {
        const b = BNDM.init("dhjalk");
        const occ = try collect(testing.allocator, &b, "dhjalkjwqnnnannanaflkjdklfj");
        defer testing.allocator.free(occ);
        try testing.expectEqualSlices(usize, &.{0}, occ);
    }
}

test "bndm max length pattern (63 symbols)" {
    var pattern: [63]u8 = undefined;
    @memset(&pattern, 'A');
    var text: [73]u8 = undefined;
    @memset(text[0..10], 'C');
    @memcpy(text[10..], &pattern);
    const b = BNDM.init(&pattern);
    const occ = try collect(testing.allocator, &b, &text);
    defer testing.allocator.free(occ);
    try testing.expectEqualSlices(usize, &.{10}, occ);
}
