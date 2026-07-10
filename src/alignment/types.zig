//! Shared alignment result types — the Zig equivalent of rust-bio's `bio_types::
//! alignment` (an external crate, so reconstructed here from its use in the
//! pairwise aligner).

const std = @import("std");

/// One step of an alignment. `xclip`/`yclip` carry the number of clipped
/// characters, so this is a tagged union rather than a plain enum.
pub const AlignmentOperation = union(enum) {
    match,
    subst,
    del,
    ins,
    xclip: usize,
    yclip: usize,
};

pub const AlignmentMode = enum { local, semiglobal, global, custom };

/// A pairwise alignment result. `operations` is heap-owned; call `deinit`.
pub const Alignment = struct {
    score: i32,
    xstart: usize,
    ystart: usize,
    xend: usize,
    yend: usize,
    xlen: usize,
    ylen: usize,
    operations: []AlignmentOperation,
    mode: AlignmentMode,

    pub fn deinit(self: *Alignment, gpa: std.mem.Allocator) void {
        gpa.free(self.operations);
        self.operations = &.{};
    }

    /// Drop all `xclip`/`yclip` operations in place (used by local/semiglobal,
    /// which report clips only implicitly via the start/end coordinates).
    pub fn filterClipOperations(self: *Alignment, gpa: std.mem.Allocator) !void {
        var w: usize = 0;
        for (self.operations) |op| {
            switch (op) {
                .match, .subst, .ins, .del => {
                    self.operations[w] = op;
                    w += 1;
                },
                .xclip, .yclip => {},
            }
        }
        self.operations = try gpa.realloc(self.operations, w);
    }
};
