const std = @import("std");
pub fn build(b: *std.Build) void {
    const target = b.standardTargetOptions(.{});
    const optimize = b.standardOptimizeOption(.{});

    // Public library module (stays dependency-free).
    const mod = b.addModule("zbio", .{
        .root_source_file = b.path("src/root.zig"),
        .target = target,
        .optimize = optimize,
    });

    const mod_tests = b.addTest(.{
        .root_module = mod,
    });
    const run_mod_tests = b.addRunArtifact(mod_tests);
    const test_step = b.step("test", "Run tests");
    test_step.dependOn(&run_mod_tests.step);
    b.installArtifact(mod_tests);

    // Benchmarks: always ReleaseFast, zBench pulled in lazily so `zig build test`
    // never compiles or links it. The bench exe imports a *separate* ReleaseFast
    // instance of zbio so the benchmarked library code is optimized (the public
    // `mod` above is Debug by default).
    const bench_step = b.step("bench", "Run zBench benchmarks (ReleaseFast)");
    const bench_optimize: std.builtin.OptimizeMode = .ReleaseFast;
    if (b.lazyDependency("zbench", .{
        .target = target,
        .optimize = bench_optimize,
    })) |zbench_dep| {
        const zbio_for_bench = b.createModule(.{
            .root_source_file = b.path("src/root.zig"),
            .target = target,
            .optimize = bench_optimize,
        });

        // Corpora rooted in src/test/ so @embedFile paths stay inside the module.
        const fixtures = b.createModule(.{
            .root_source_file = b.path("src/test/fixtures.zig"),
            .target = target,
            .optimize = bench_optimize,
        });

        const bench_mod = b.createModule(.{
            .root_source_file = b.path("bench/main.zig"),
            .target = target,
            .optimize = bench_optimize,
            .imports = &.{
                .{ .name = "zbio", .module = zbio_for_bench },
                .{ .name = "zbench", .module = zbench_dep.module("zbench") },
                .{ .name = "fixtures", .module = fixtures },
            },
        });

        const bench_exe = b.addExecutable(.{
            .name = "zbio-bench",
            .root_module = bench_mod,
        });

        const run_bench = b.addRunArtifact(bench_exe);
        bench_step.dependOn(&run_bench.step);
    }

    // Head-to-head comparison harness (vs rust-bio). Always ReleaseFast and links
    // libc so the Zig side hits the system malloc, matching Rust's allocator. The
    // binary is installed; bench/compare/run.sh invokes it with a corpus-path arg.
    const compare_step = b.step("compare", "Build the rust-bio comparison harness (ReleaseFast)");
    const compare_optimize: std.builtin.OptimizeMode = .ReleaseFast;
    const zbio_for_compare = b.createModule(.{
        .root_source_file = b.path("src/root.zig"),
        .target = target,
        .optimize = compare_optimize,
    });
    const compare_mod = b.createModule(.{
        .root_source_file = b.path("bench/compare/zig/main.zig"),
        .target = target,
        .optimize = compare_optimize,
        .link_libc = true,
        .imports = &.{
            .{ .name = "zbio", .module = zbio_for_compare },
        },
    });
    const compare_exe = b.addExecutable(.{
        .name = "zbio-compare",
        .root_module = compare_mod,
    });
    const install_compare = b.addInstallArtifact(compare_exe, .{});
    compare_step.dependOn(&install_compare.step);
}
