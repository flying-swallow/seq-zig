//! Shift-And bit-parallel exact matching. Ported from rust-bio's
//! `pattern_matching::shift_and`. Patterns must have fewer than 64 symbols
//! (state set packed in a u64). Complexity O(n). No allocation.

const std = @import("std");
const testing = std.testing;

pub const Masks = struct {
    table: [256]u64,
    accept: u64,
};

/// Build the per-symbol bit masks and accept mask. `reverse` iterates the
/// pattern back-to-front (BNDM uses this over the reversed pattern).
pub fn buildMasks(pattern: []const u8, comptime reverse: bool) Masks {
    var m = Masks{ .table = @splat(0), .accept = 0 };
    var bit: u64 = 1;
    if (reverse) {
        var idx = pattern.len;
        while (idx > 0) {
            idx -= 1;
            m.table[pattern[idx]] |= bit;
            bit <<= 1;
        }
    } else {
        for (pattern) |c| {
            m.table[c] |= bit;
            bit <<= 1;
        }
    }
    m.accept = bit >> 1; // bit for the last pattern position
    return m;
}

/// Convenience: forward masks.
pub fn masks(pattern: []const u8) Masks {
    return buildMasks(pattern, false);
}

pub const ShiftAnd = struct {
    m: usize,
    table: [256]u64,
    accept: u64,

    pub fn init(pattern: []const u8) ShiftAnd {
        std.debug.assert(pattern.len < 64); // pattern of less than 64 symbols
        const mk = masks(pattern);
        return .{ .m = pattern.len, .table = mk.table, .accept = mk.accept };
    }

    pub fn findAll(self: *const ShiftAnd, text: []const u8) Matches {
        return .{ .shiftand = self, .active = 0, .text = text, .pos = 0 };
    }

    pub const Matches = struct {
        shiftand: *const ShiftAnd,
        active: u64,
        text: []const u8,
        pos: usize,

        pub fn next(self: *Matches) ?usize {
            const sa = self.shiftand;
            while (self.pos < self.text.len) {
                const i = self.pos;
                const c = self.text[i];
                self.pos += 1;
                self.active = ((self.active << 1) | 1) & sa.table[c];
                if (self.active & sa.accept != 0) return i + 1 - sa.m;
            }
            return null;
        }
    };
};

// Tests -----------------------------------------------------------------------

fn collect(gpa: std.mem.Allocator, sa: *const ShiftAnd, text: []const u8) ![]usize {
    var out: std.ArrayListUnmanaged(usize) = .empty;
    var it = sa.findAll(text);
    while (it.next()) |p| try out.append(gpa, p);
    return out.toOwnedSlice(gpa);
}

test "shift_and find_all" {
    const sa = ShiftAnd.init("qnnnannan");
    const occ = try collect(testing.allocator, &sa, "dhjalkjwqnnnannanaflkjdklfj");
    defer testing.allocator.free(occ);
    try testing.expectEqualSlices(usize, &.{8}, occ);
}

test "shift_and issue_416 (match at 0)" {
    const sa = ShiftAnd.init("CC");
    const occ = try collect(testing.allocator, &sa, "CCTTTTTTTTTTTTTTT");
    defer testing.allocator.free(occ);
    try testing.expectEqualSlices(usize, &.{0}, occ);
}

test "shift_and multiple finds" {
    const sa = ShiftAnd.init("CC");
    const occ = try collect(testing.allocator, &sa, "CCTCCTCC");
    defer testing.allocator.free(occ);
    try testing.expectEqualSlices(usize, &.{ 0, 3, 6 }, occ);
}

test "shift_and max length pattern (63 symbols)" {
    var pattern: [63]u8 = undefined;
    @memset(&pattern, 'A');
    var text: [73]u8 = undefined;
    @memset(text[0..10], 'C');
    @memcpy(text[10..], &pattern);
    const sa = ShiftAnd.init(&pattern);
    const occ = try collect(testing.allocator, &sa, &text);
    defer testing.allocator.free(occ);
    try testing.expectEqualSlices(usize, &.{10}, occ);
}
