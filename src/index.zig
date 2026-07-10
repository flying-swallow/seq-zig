//! Full-text index namespace: suffix array, BWT, and FM-index, ported from
//! rust-bio's `data_structures`. (The bidirectional FMD-index and sampled
//! suffix array are deferred.)

pub const suffix_array = @import("index/suffix_array.zig");
pub const bwt = @import("index/bwt.zig");
pub const fmindex = @import("index/fmindex.zig");

test {
    _ = suffix_array;
    _ = bwt;
    _ = fmindex;
}
