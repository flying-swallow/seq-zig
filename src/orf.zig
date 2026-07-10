//! One-way (forward-strand) open reading frame finder, ported from rust-bio's
//! `seq_analysis::orf`. A streaming state machine keeps, per reading frame, a
//! stack of open start positions; on a stop codon it emits every sufficiently
//! long ORF for that frame (in ascending start order) into a FIFO. To scan the
//! reverse strand, reverse-complement the sequence and run again.

const std = @import("std");
const testing = std.testing;

/// An ORF: half-open byte range [start, end) and reading-frame offset (0,1,2).
pub const Orf = struct {
    start: usize,
    end: usize,
    offset: i8,
};

pub const Finder = struct {
    start_codons: []const [3]u8,
    stop_codons: []const [3]u8,
    min_len: usize,

    pub fn init(start_codons: []const [3]u8, stop_codons: []const [3]u8, min_len: usize) Finder {
        return .{ .start_codons = start_codons, .stop_codons = stop_codons, .min_len = min_len };
    }

    /// Streaming iterator over ORFs in `seq`. Caller must `deinit` it.
    pub fn findAll(self: *const Finder, gpa: std.mem.Allocator, seq: []const u8) Matches {
        return .{
            .finder = self,
            .gpa = gpa,
            .seq = seq,
            .i = 0,
            .codon = undefined,
            .clen = 0,
            .start_pos = .{ .empty, .empty, .empty },
            .found = .empty,
            .found_head = 0,
        };
    }

    /// Collect every ORF into an owned slice (convenience for callers/tests).
    pub fn collectAll(self: *const Finder, gpa: std.mem.Allocator, seq: []const u8) ![]Orf {
        var m = self.findAll(gpa, seq);
        defer m.deinit();
        var out: std.ArrayListUnmanaged(Orf) = .empty;
        errdefer out.deinit(gpa);
        while (try m.next()) |orf| try out.append(gpa, orf);
        return out.toOwnedSlice(gpa);
    }
};

fn codonMatches(list: []const [3]u8, codon: [3]u8) bool {
    for (list) |c| {
        if (std.mem.eql(u8, &c, &codon)) return true;
    }
    return false;
}

pub const Matches = struct {
    finder: *const Finder,
    gpa: std.mem.Allocator,
    seq: []const u8,
    i: usize,
    codon: [3]u8,
    clen: u8, // number of valid bytes in the sliding codon window (0..3)
    start_pos: [3]std.ArrayListUnmanaged(usize),
    found: std.ArrayListUnmanaged(Orf), // FIFO drained via found_head
    found_head: usize,

    pub fn deinit(self: *Matches) void {
        for (&self.start_pos) |*sp| sp.deinit(self.gpa);
        self.found.deinit(self.gpa);
        self.* = undefined;
    }

    fn dequeue(self: *Matches) ?Orf {
        if (self.found_head < self.found.items.len) {
            const orf = self.found.items[self.found_head];
            self.found_head += 1;
            return orf;
        }
        return null;
    }

    pub fn next(self: *Matches) !?Orf {
        // Return any ORFs already queued.
        if (self.dequeue()) |orf| return orf;

        while (self.i < self.seq.len) {
            const index = self.i;
            const nuc = self.seq[index];
            self.i += 1;

            // Slide the 3-byte codon window.
            if (self.clen == 3) {
                self.codon[0] = self.codon[1];
                self.codon[1] = self.codon[2];
                self.codon[2] = nuc;
            } else {
                self.codon[self.clen] = nuc;
                self.clen += 1;
            }
            const have_codon = self.clen == 3;
            const offset = (index + 1) % 3;

            // Entering an ORF?
            if (have_codon and codonMatches(self.finder.start_codons, self.codon)) {
                try self.start_pos[offset].append(self.gpa, index);
            }
            // Inside an ORF for this frame — check for a stop.
            if (self.start_pos[offset].items.len > 0 and
                have_codon and codonMatches(self.finder.stop_codons, self.codon))
            {
                for (self.start_pos[offset].items) |sp| {
                    if (index + 1 - sp > self.finder.min_len) {
                        try self.found.append(self.gpa, .{
                            .start = sp - 2,
                            .end = index + 1,
                            .offset = @intCast(offset),
                        });
                    } else {
                        // starts are ascending; the rest are shorter still.
                        break;
                    }
                }
                self.start_pos[offset].clearRetainingCapacity();
            }

            if (self.dequeue()) |orf| return orf;
        }
        return null;
    }
};

// Tests -----------------------------------------------------------------------

const test_starts = [_][3]u8{"ATG".*};
const test_stops = [_][3]u8{ "TGA".*, "TAG".*, "TAA".* };

fn basicFinder() Finder {
    return Finder.init(&test_starts, &test_stops, 5);
}

test "orf: no orf" {
    const f = basicFinder();
    const orfs = try f.collectAll(testing.allocator, "ACGGCTAGAAAAGGCTAGAAAA");
    defer testing.allocator.free(orfs);
    try testing.expectEqual(@as(usize, 0), orfs.len);
}

test "orf: one orf no offset" {
    const f = basicFinder();
    const orfs = try f.collectAll(testing.allocator, "GGGATGGGGTGAGGG");
    defer testing.allocator.free(orfs);
    try testing.expectEqualSlices(Orf, &.{.{ .start = 3, .end = 12, .offset = 0 }}, orfs);
}

test "orf: one orf with offset" {
    const f = basicFinder();
    const orfs = try f.collectAll(testing.allocator, "AGGGATGGGGTGAGGG");
    defer testing.allocator.free(orfs);
    try testing.expectEqualSlices(Orf, &.{.{ .start = 4, .end = 13, .offset = 1 }}, orfs);
}

test "orf: two orfs different offsets" {
    const f = basicFinder();
    const orfs = try f.collectAll(testing.allocator, "ATGGGGTGAGGGGGATGGAAAAATAAG");
    defer testing.allocator.free(orfs);
    try testing.expectEqualSlices(Orf, &.{
        .{ .start = 0, .end = 9, .offset = 0 },
        .{ .start = 14, .end = 26, .offset = 2 },
    }, orfs);
}

test "orf: three nested and offset orfs" {
    const f = basicFinder();
    const orfs = try f.collectAll(testing.allocator, "ATGGGGATGGGGGGATGGAAAAATAAGTAG");
    defer testing.allocator.free(orfs);
    try testing.expectEqualSlices(Orf, &.{
        .{ .start = 14, .end = 26, .offset = 2 },
        .{ .start = 0, .end = 30, .offset = 0 },
        .{ .start = 6, .end = 30, .offset = 0 },
    }, orfs);
}
