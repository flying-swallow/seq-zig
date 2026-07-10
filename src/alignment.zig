//! Alignment namespace aggregator. `align` is a reserved Zig keyword, so this
//! module is named `alignment` (matching rust-bio's `src/alignment/`).

pub const types = @import("alignment/types.zig");
pub const distance = @import("alignment/distance.zig");
pub const pairwise = @import("alignment/pairwise.zig");

// Re-export the shared alignment result types at the namespace root.
pub const Alignment = types.Alignment;
pub const AlignmentOperation = types.AlignmentOperation;
pub const AlignmentMode = types.AlignmentMode;

test {
    _ = types;
    _ = distance;
    _ = pairwise;
}
