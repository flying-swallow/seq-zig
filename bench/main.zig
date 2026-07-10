//! zbio benchmark suite entry point. Built only by `zig build bench` (ReleaseFast).
//!
//! Each domain file exposes `pub fn register(bench: *zbench.Benchmark) !void`;
//! adding a benchmark is a one-line `bench.add(...)` in the relevant file, no
//! build.zig change required.

const std = @import("std");
const zbench = @import("zbench");

const io_bench = @import("io.zig");
const algo_bench = @import("algorithms.zig");

pub fn main(init: std.process.Init) !void {
    const io = init.io;
    const stdout: std.Io.File = .stdout();

    var bench = zbench.Benchmark.init(init.gpa, .{});
    defer bench.deinit();

    try io_bench.register(&bench);
    try algo_bench.register(&bench);

    try bench.run(io, stdout);
}
