//! Backward Oracle Matching (BOM). Ported from rust-bio's
//! `pattern_matching::bom`. rust uses `Vec<VecMap<usize>>` for the oracle's
//! transitions; here it's a dense `[]usize` of `m * (maxsym+1)` with 0 meaning
//! "no edge" (all real edge targets are >= 1). Best case O(n/m). Owns the
//! transition table — call `deinit`.

const std = @import("std");
const testing = std.testing;

pub const BOM = struct {
    m: usize,
    maxsym: usize,
    stride: usize,
    table: []usize, // m rows x stride cols; 0 = no edge

    pub fn init(gpa: std.mem.Allocator, pattern: []const u8) !BOM {
        const m = pattern.len;
        std.debug.assert(m > 0); // non-empty pattern
        var maxsym: usize = 0;
        for (pattern) |c| maxsym = @max(maxsym, c);
        const stride = maxsym + 1;

        const table = try gpa.alloc(usize, m * stride);
        errdefer gpa.free(table);
        @memset(table, 0);

        // suff[i] = state where the longest suffix of pattern[..i+1] not ending
        // in i ends; None initially.
        const suff = try gpa.alloc(?usize, m + 1);
        defer gpa.free(suff);
        @memset(suff, null);

        var j: usize = 0;
        while (j < m) : (j += 1) {
            const a: usize = pattern[m - 1 - j]; // reversed pattern
            const i = j + 1;
            table[j * stride + a] = i; // inner edge: reading a leads to state i
            var k = suff[i - 1];
            while (k) |k_| {
                if (table[k_ * stride + a] != 0) break;
                table[k_ * stride + a] = i;
                k = suff[k_];
            }
            suff[i] = if (k) |k_| table[k_ * stride + a] else 0;
        }

        return .{ .m = m, .maxsym = maxsym, .stride = stride, .table = table };
    }

    pub fn deinit(self: *BOM, gpa: std.mem.Allocator) void {
        gpa.free(self.table);
        self.* = undefined;
    }

    pub fn delta(self: *const BOM, q: usize, a: u8) ?usize {
        if (q >= self.m) return null;
        if (a > self.maxsym) return null;
        const v = self.table[q * self.stride + a];
        return if (v == 0) null else v;
    }

    pub fn findAll(self: *const BOM, text: []const u8) Matches {
        return .{ .bom = self, .text = text, .window = self.m };
    }

    pub const Matches = struct {
        bom: *const BOM,
        text: []const u8,
        window: usize,

        pub fn next(self: *Matches) ?usize {
            const b = self.bom;
            while (self.window <= self.text.len) {
                var q: ?usize = 0;
                var j: usize = 1;
                while (j <= b.m) {
                    if (q) |q_| {
                        q = b.delta(q_, self.text[self.window - j]);
                        j += 1;
                    } else break;
                }
                const i = self.window - b.m; // putative start
                self.window += b.m + 2 - j;
                if (q != null) return i;
            }
            return null;
        }
    };
};

// Tests -----------------------------------------------------------------------

fn collect(gpa: std.mem.Allocator, b: *const BOM, text: []const u8) ![]usize {
    var out: std.ArrayListUnmanaged(usize) = .empty;
    var it = b.findAll(text);
    while (it.next()) |p| try out.append(gpa, p);
    return out.toOwnedSlice(gpa);
}

test "bom delta" {
    var b = try BOM.init(testing.allocator, "qnnnannan");
    defer b.deinit(testing.allocator);
    try testing.expectEqual(@as(?usize, 1), b.delta(0, 'n'));
    try testing.expectEqual(@as(?usize, 2), b.delta(1, 'a'));
    try testing.expectEqual(@as(?usize, 3), b.delta(2, 'n'));
    try testing.expectEqual(@as(?usize, 4), b.delta(3, 'n'));
    try testing.expectEqual(@as(?usize, 5), b.delta(4, 'a'));
    try testing.expectEqual(@as(?usize, 6), b.delta(5, 'n'));
    try testing.expectEqual(@as(?usize, 7), b.delta(6, 'n'));
    try testing.expectEqual(@as(?usize, 8), b.delta(7, 'n'));
    try testing.expectEqual(@as(?usize, 9), b.delta(8, 'q'));
    try testing.expectEqual(@as(?usize, 2), b.delta(0, 'a'));
    try testing.expectEqual(@as(?usize, 9), b.delta(0, 'q'));
    try testing.expectEqual(@as(?usize, 4), b.delta(1, 'n'));
    try testing.expectEqual(@as(?usize, 9), b.delta(1, 'q'));
    try testing.expectEqual(@as(?usize, 8), b.delta(4, 'n'));
    try testing.expectEqual(@as(?usize, 9), b.delta(4, 'q'));
    try testing.expectEqual(@as(?usize, null), b.delta(9, 'a')); // q >= m
}

test "bom find_all" {
    {
        var b = try BOM.init(testing.allocator, "GAAAA");
        defer b.deinit(testing.allocator);
        const occ = try collect(testing.allocator, &b, "ACGGCTAGGAAAAAGACTGAGGACTGAAAA");
        defer testing.allocator.free(occ);
        try testing.expectEqualSlices(usize, &.{ 8, 25 }, occ);
    }
    {
        var b = try BOM.init(testing.allocator, "qnnnannan");
        defer b.deinit(testing.allocator);
        const occ = try collect(testing.allocator, &b, "dhjalkjwqnnnannanaflkjdklfj");
        defer testing.allocator.free(occ);
        try testing.expectEqualSlices(usize, &.{8}, occ);
    }
    {
        var b = try BOM.init(testing.allocator, "dhjalk");
        defer b.deinit(testing.allocator);
        const occ = try collect(testing.allocator, &b, "dhjalkjwqnnnannanaflkjdklfj");
        defer testing.allocator.free(occ);
        try testing.expectEqualSlices(usize, &.{0}, occ);
    }
}
