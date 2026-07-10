//! Exact pattern-matching algorithms, ported from rust-bio's `pattern_matching`.
//! Each matcher exposes `init` + `findAll(text)` returning a `Matches` iterator
//! whose `next()` yields start positions.

pub const horspool = @import("pm/horspool.zig");
pub const kmp = @import("pm/kmp.zig");
pub const shift_and = @import("pm/shift_and.zig");
pub const bndm = @import("pm/bndm.zig");
pub const bom = @import("pm/bom.zig");

test {
    _ = horspool;
    _ = kmp;
    _ = shift_and;
    _ = bndm;
    _ = bom;
}
