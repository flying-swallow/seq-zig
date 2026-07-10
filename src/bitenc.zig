//! Fixed-width bit encoding, ported from rust-bio's `data_structures::bitenc`.
//! Packs a sequence of small values (1..=8 bits each) into u32 blocks. Pairs
//! with `alphabet.RankTransform` to store rank-encoded sequences (e.g. 2-bit
//! DNA). `usable_bits_per_block = 32 - 32 % width`, so for widths that don't
//! divide 32 some high bits of each block go unused.

const std = @import("std");
const testing = std.testing;

fn widthMask(width: u5) u32 {
    return (@as(u32, 1) << width) - 1;
}

pub const BitEnc = struct {
    storage: std.ArrayListUnmanaged(u32) = .empty,
    width: u5, // 1..=8
    mask: u32,
    len: usize,
    usable_bits_per_block: u6, // up to 32

    /// Create an encoder with the given bit width (1..=8).
    pub fn init(width: u5) BitEnc {
        std.debug.assert(width >= 1 and width <= 8);
        return .{
            .storage = .empty,
            .width = width,
            .mask = widthMask(width),
            .len = 0,
            .usable_bits_per_block = @intCast(32 - 32 % @as(usize, width)),
        };
    }

    pub fn deinit(self: *BitEnc, gpa: std.mem.Allocator) void {
        self.storage.deinit(gpa);
        self.* = undefined;
    }

    // Bit positions within a block are 0..31, but the pushValues fill loop can
    // momentarily step a counter past 31, so bit positions are carried as usize
    // and narrowed to u5 (the shift-amount type for u32) at each shift site.
    fn addr(self: *const BitEnc, i: usize) struct { block: usize, bit: usize } {
        const k = i * self.width;
        return .{
            .block = k / self.usable_bits_per_block,
            .bit = k % self.usable_bits_per_block,
        };
    }

    fn getByAddr(self: *const BitEnc, block: usize, bit: usize) u8 {
        return @intCast((self.storage.items[block] >> @intCast(bit)) & self.mask);
    }

    fn setByAddr(self: *BitEnc, block: usize, bit: usize, value: u8) void {
        const sh: u5 = @intCast(bit);
        const m = self.mask << sh;
        self.storage.items[block] |= m;
        self.storage.items[block] ^= m;
        self.storage.items[block] |= ((@as(u32, value) & self.mask) << sh);
    }

    /// Append one value. Complexity O(1) amortized.
    pub fn push(self: *BitEnc, gpa: std.mem.Allocator, value: u8) !void {
        const a = self.addr(self.len);
        if (a.bit == 0) try self.storage.append(gpa, 0);
        self.setByAddr(a.block, a.bit, value);
        self.len += 1;
    }

    /// Append `value` `n` times. Complexity O(n).
    pub fn pushValues(self: *BitEnc, gpa: std.mem.Allocator, n_in: usize, value: u8) !void {
        var n = n_in;

        // Fill up any remaining free slots in the current (partial) block.
        {
            const a = self.addr(self.len);
            if (a.bit > 0) {
                var bit: usize = a.bit;
                while (bit < 32 and n > 0) : (bit += self.width) {
                    self.setByAddr(a.block, bit, value);
                    n -= 1;
                    self.len += 1;
                }
            }
        }

        if (n > 0) {
            // Build a block completely filled with copies of `value`.
            var value_block: u32 = 0;
            {
                var v: u32 = value;
                var c: usize = 0;
                while (c < 32 / @as(usize, self.width)) : (c += 1) {
                    value_block |= v;
                    v <<= self.width;
                }
            }

            // Grow storage to `a.block` full value-blocks, then a partial block.
            // Vec::resize(block, value_block) fills newly added slots; ArrayList
            // resize leaves them undefined, so memset the grown region.
            const i = self.len + n;
            const a = self.addr(i);
            const old_len = self.storage.items.len;
            try self.storage.resize(gpa, a.block);
            if (a.block > old_len) @memset(self.storage.items[old_len..a.block], value_block);

            if (a.bit > 0) {
                const shift_amt: u5 = @intCast(@as(usize, self.usable_bits_per_block) - a.bit);
                try self.storage.append(gpa, value_block >> shift_amt);
            }
            self.len = i;
        }
    }

    /// Replace the value at position `i`.
    pub fn set(self: *BitEnc, i: usize, value: u8) void {
        const a = self.addr(i);
        self.setByAddr(a.block, a.bit, value);
    }

    /// Get the value at position `i`, or null if out of range.
    pub fn get(self: *const BitEnc, i: usize) ?u8 {
        if (i >= self.len) return null;
        const a = self.addr(i);
        return self.getByAddr(a.block, a.bit);
    }

    pub fn clear(self: *BitEnc) void {
        self.storage.clearRetainingCapacity();
        self.len = 0;
    }

    pub fn nrBlocks(self: *const BitEnc) usize {
        return self.storage.items.len;
    }
    pub fn nrSymbols(self: *const BitEnc) usize {
        return self.len;
    }
    pub fn isEmpty(self: *const BitEnc) bool {
        return self.len == 0;
    }

    pub fn iterator(self: *const BitEnc) Iterator {
        return .{ .bitenc = self, .i = 0 };
    }

    pub const Iterator = struct {
        bitenc: *const BitEnc,
        i: usize,
        pub fn next(self: *Iterator) ?u8 {
            const v = self.bitenc.get(self.i);
            self.i += 1;
            return v;
        }
    };
};

// Tests -----------------------------------------------------------------------

fn collect(gpa: std.mem.Allocator, be: *const BitEnc) ![]u8 {
    var out: std.ArrayListUnmanaged(u8) = .empty;
    var it = be.iterator();
    while (it.next()) |v| try out.append(gpa, v);
    return out.toOwnedSlice(gpa);
}

test "bitenc push/set/get" {
    var be = BitEnc.init(2);
    defer be.deinit(testing.allocator);
    try be.push(testing.allocator, 0);
    try be.push(testing.allocator, 2);
    try be.push(testing.allocator, 1);
    {
        const vals = try collect(testing.allocator, &be);
        defer testing.allocator.free(vals);
        try testing.expectEqualSlices(u8, &.{ 0, 2, 1 }, vals);
    }
    be.set(1, 3);
    {
        const vals = try collect(testing.allocator, &be);
        defer testing.allocator.free(vals);
        try testing.expectEqualSlices(u8, &.{ 0, 3, 1 }, vals);
    }
}

test "bitenc pushValues edge cases (width 7)" {
    var be = BitEnc.init(7);
    defer be.deinit(testing.allocator);
    // width 7 -> 4 values per block, 4 leftover bits.
    try be.pushValues(testing.allocator, 5, 42);
    {
        const vals = try collect(testing.allocator, &be);
        defer testing.allocator.free(vals);
        try testing.expectEqualSlices(u8, &.{ 42, 42, 42, 42, 42 }, vals);
    }
    try testing.expectEqual(@as(usize, 2), be.nrBlocks());
    try testing.expectEqual(@as(usize, 5), be.nrSymbols());

    try be.pushValues(testing.allocator, 1, 23);
    {
        const vals = try collect(testing.allocator, &be);
        defer testing.allocator.free(vals);
        try testing.expectEqualSlices(u8, &.{ 42, 42, 42, 42, 42, 23 }, vals);
    }

    try be.pushValues(testing.allocator, 12, 17);
    {
        const vals = try collect(testing.allocator, &be);
        defer testing.allocator.free(vals);
        try testing.expectEqualSlices(u8, &.{
            42, 42, 42, 42, 42, 23, 17, 17, 17,
            17, 17, 17, 17, 17, 17, 17, 17, 17,
        }, vals);
    }
    try testing.expectEqual(@as(usize, 5), be.nrBlocks());
    try testing.expectEqual(@as(usize, 18), be.nrSymbols());
}

test "bitenc pushValues 32 zeros (width 2)" {
    var be = BitEnc.init(2);
    defer be.deinit(testing.allocator);
    try be.pushValues(testing.allocator, 32, 0);
    try testing.expectEqual(@as(usize, 2), be.nrBlocks());
    try testing.expectEqual(@as(usize, 32), be.nrSymbols());
}

test "bitenc width 2..8 stress (issue29)" {
    var w: u5 = 2;
    while (w < 9) : (w += 1) {
        var be = BitEnc.init(w);
        defer be.deinit(testing.allocator);
        var i: usize = 0;
        while (i < 1000) : (i += 1) try be.push(testing.allocator, 1);
        try testing.expectEqual(@as(usize, 1000), be.nrSymbols());
        try testing.expectEqual(@as(?u8, 1), be.get(999));
    }
}
