//! Pairwise alignment via a generalized Smith-Waterman (Gotoh affine-gap)
//! algorithm, ported from rust-bio's `alignment::pairwise`. A single `custom`
//! core parameterized by four clip penalties yields global / semiglobal / local
//! alignments; O(n*m) time, O(n) score-column space (plus the full traceback
//! matrix). The scoring function is a comptime type exposing `score(a, b) i32`
//! (rust's `MatchFunc` trait) — use `MatchParams` for constant match/mismatch
//! or `FnMatch(f)` to wrap a matrix function such as `scores.blosum62`.

const std = @import("std");
const testing = std.testing;
const types = @import("types.zig");

pub const Alignment = types.Alignment;
pub const AlignmentOperation = types.AlignmentOperation;
pub const AlignmentMode = types.AlignmentMode;

/// 'Negative infinity' score, ~0.4 * i32.MIN. Chosen so that sums like
/// `MIN_SCORE + MIN_SCORE` or `MIN_SCORE + gap_extend*(k-1)` don't underflow
/// i32 (which panics in Zig safe builds).
pub const MIN_SCORE: i32 = -858_993_459;

// Traceback move codes (must fit in 4 bits; TB_MAX guards set_bits).
const TB_START: u16 = 0b0000;
const TB_INS: u16 = 0b0001;
const TB_DEL: u16 = 0b0010;
const TB_SUBST: u16 = 0b0011;
const TB_MATCH: u16 = 0b0100;
const TB_XCLIP_PREFIX: u16 = 0b0101;
const TB_XCLIP_SUFFIX: u16 = 0b0110;
const TB_YCLIP_PREFIX: u16 = 0b0111;
const TB_YCLIP_SUFFIX: u16 = 0b1000;
const TB_MAX: u16 = 0b1000;

const I_POS: u4 = 0;
const D_POS: u4 = 4;
const S_POS: u4 = 8;

/// Packs the I, D and S traceback moves (4 bits each) into one u16.
const TracebackCell = struct {
    v: u16 = 0,

    fn setBits(self: *TracebackCell, pos: u4, value: u16) void {
        std.debug.assert(value <= TB_MAX);
        const bits: u16 = @as(u16, 0b1111) << pos;
        self.v = (self.v & ~bits) | (value << pos);
    }
    fn setIBits(self: *TracebackCell, value: u16) void {
        self.setBits(I_POS, value);
    }
    fn setDBits(self: *TracebackCell, value: u16) void {
        self.setBits(D_POS, value);
    }
    fn setSBits(self: *TracebackCell, value: u16) void {
        self.setBits(S_POS, value);
    }
    fn getBits(self: TracebackCell, pos: u4) u16 {
        return (self.v >> pos) & 0b1111;
    }
    fn getIBits(self: TracebackCell) u16 {
        return self.getBits(I_POS);
    }
    fn getDBits(self: TracebackCell) u16 {
        return self.getBits(D_POS);
    }
    fn getSBits(self: TracebackCell) u16 {
        return self.getBits(S_POS);
    }
    fn setAll(self: *TracebackCell, value: u16) void {
        self.setIBits(value);
        self.setDBits(value);
        self.setSBits(value);
    }
};

const Traceback = struct {
    rows: usize = 0,
    cols: usize = 0,
    matrix: std.ArrayListUnmanaged(TracebackCell) = .empty,

    fn deinit(self: *Traceback, gpa: std.mem.Allocator) void {
        self.matrix.deinit(gpa);
    }
    fn init(self: *Traceback, gpa: std.mem.Allocator, m: usize, n: usize) !void {
        self.rows = m + 1;
        self.cols = n + 1;
        try self.matrix.resize(gpa, self.rows * self.cols);
        @memset(self.matrix.items, .{ .v = 0 }); // all cells TB_START (0)
    }
    fn set(self: *Traceback, i: usize, j: usize, v: TracebackCell) void {
        self.matrix.items[i * self.cols + j] = v;
    }
    fn get(self: *const Traceback, i: usize, j: usize) TracebackCell {
        return self.matrix.items[i * self.cols + j];
    }
    fn getPtr(self: *Traceback, i: usize, j: usize) *TracebackCell {
        return &self.matrix.items[i * self.cols + j];
    }
};

/// Constant match/mismatch scoring (rust's `MatchParams`).
pub const MatchParams = struct {
    match_score: i32,
    mismatch_score: i32,

    pub fn init(match_score: i32, mismatch_score: i32) MatchParams {
        std.debug.assert(match_score >= 0);
        std.debug.assert(mismatch_score <= 0);
        return .{ .match_score = match_score, .mismatch_score = mismatch_score };
    }
    pub fn score(self: MatchParams, a: u8, b: u8) i32 {
        return if (a == b) self.match_score else self.mismatch_score;
    }
};

/// Wrap a plain `fn(u8, u8) i32` (e.g. `scores.blosum62`) as a MatchFunc type.
pub fn FnMatch(comptime f: fn (u8, u8) i32) type {
    return struct {
        pub fn score(_: @This(), a: u8, b: u8) i32 {
            return f(a, b);
        }
    };
}

/// Affine-gap scoring parameters, generic over the match-function type.
pub fn Scoring(comptime MatchFunc: type) type {
    return struct {
        const Self = @This();

        gap_open: i32,
        gap_extend: i32,
        match_fn: MatchFunc,
        match_scores: ?[2]i32 = null,
        xclip_prefix: i32 = MIN_SCORE,
        xclip_suffix: i32 = MIN_SCORE,
        yclip_prefix: i32 = MIN_SCORE,
        yclip_suffix: i32 = MIN_SCORE,

        /// New scoring with clip penalties defaulting to MIN_SCORE (i.e. global).
        pub fn init(gap_open: i32, gap_extend: i32, match_fn: MatchFunc) Self {
            std.debug.assert(gap_open <= 0);
            std.debug.assert(gap_extend <= 0);
            return .{ .gap_open = gap_open, .gap_extend = gap_extend, .match_fn = match_fn };
        }

        pub fn xclip(self: Self, penalty: i32) Self {
            var s = self;
            s.xclip_prefix = penalty;
            s.xclip_suffix = penalty;
            return s;
        }
        pub fn xclipPrefix(self: Self, penalty: i32) Self {
            var s = self;
            s.xclip_prefix = penalty;
            return s;
        }
        pub fn xclipSuffix(self: Self, penalty: i32) Self {
            var s = self;
            s.xclip_suffix = penalty;
            return s;
        }
        pub fn yclip(self: Self, penalty: i32) Self {
            var s = self;
            s.yclip_prefix = penalty;
            s.yclip_suffix = penalty;
            return s;
        }
        pub fn yclipPrefix(self: Self, penalty: i32) Self {
            var s = self;
            s.yclip_prefix = penalty;
            return s;
        }
        pub fn yclipSuffix(self: Self, penalty: i32) Self {
            var s = self;
            s.yclip_suffix = penalty;
            return s;
        }
    };
}

/// Build constant-score scoring (rust's `Scoring::from_scores`), clips MIN_SCORE.
pub fn scoringFromScores(gap_open: i32, gap_extend: i32, match_score: i32, mismatch_score: i32) Scoring(MatchParams) {
    var s = Scoring(MatchParams).init(gap_open, gap_extend, MatchParams.init(match_score, mismatch_score));
    s.match_scores = .{ match_score, mismatch_score };
    return s;
}

pub fn Aligner(comptime MatchFunc: type) type {
    return struct {
        const Self = @This();
        const Score = std.ArrayListUnmanaged(i32);

        gpa: std.mem.Allocator,
        I: [2]Score = .{ .empty, .empty },
        D: [2]Score = .{ .empty, .empty },
        S: [2]Score = .{ .empty, .empty },
        Lx: std.ArrayListUnmanaged(usize) = .empty,
        Ly: std.ArrayListUnmanaged(usize) = .empty,
        Sn: Score = .empty,
        traceback: Traceback = .{},
        scoring: Scoring(MatchFunc),

        pub fn init(gpa: std.mem.Allocator, gap_open: i32, gap_extend: i32, match_fn: MatchFunc) Self {
            return .{ .gpa = gpa, .scoring = Scoring(MatchFunc).init(gap_open, gap_extend, match_fn) };
        }

        pub fn withScoring(gpa: std.mem.Allocator, scoring: Scoring(MatchFunc)) Self {
            return .{ .gpa = gpa, .scoring = scoring };
        }

        pub fn deinit(self: *Self) void {
            for (&self.I) |*l| l.deinit(self.gpa);
            for (&self.D) |*l| l.deinit(self.gpa);
            for (&self.S) |*l| l.deinit(self.gpa);
            self.Lx.deinit(self.gpa);
            self.Ly.deinit(self.gpa);
            self.Sn.deinit(self.gpa);
            self.traceback.deinit(self.gpa);
            self.* = undefined;
        }

        fn fillI32(self: *Self, list: *Score, len: usize, val: i32) !void {
            try list.resize(self.gpa, len);
            @memset(list.items, val);
        }
        fn fillUsize(self: *Self, list: *std.ArrayListUnmanaged(usize), len: usize, val: usize) !void {
            try list.resize(self.gpa, len);
            @memset(list.items, val);
        }

        /// The generalized SW core. Returns an owned Alignment (caller deinits).
        pub fn custom(self: *Self, x: []const u8, y: []const u8) !Alignment {
            const m = x.len;
            const n = y.len;
            const sc = &self.scoring;
            try self.traceback.init(self.gpa, m, n);

            // Initial conditions (column j=0) for both parity slots.
            var k: usize = 0;
            while (k < 2) : (k += 1) {
                try self.fillI32(&self.D[k], m + 1, MIN_SCORE);
                try self.fillI32(&self.I[k], m + 1, MIN_SCORE);
                try self.fillI32(&self.S[k], m + 1, MIN_SCORE);
                self.S[k].items[0] = 0;

                if (k == 0) {
                    var tb0 = TracebackCell{};
                    tb0.setAll(TB_START);
                    self.traceback.set(0, 0, tb0);
                    try self.fillUsize(&self.Lx, n + 1, 0);
                    try self.fillUsize(&self.Ly, m + 1, 0);
                    try self.fillI32(&self.Sn, m + 1, MIN_SCORE);
                    self.Sn.items[0] = sc.yclip_suffix;
                    self.Ly.items[0] = n;
                }

                const Ik = self.I[k].items;
                const Sk = self.S[k].items;
                var i: usize = 1;
                while (i <= m) : (i += 1) {
                    var tb = TracebackCell{};
                    tb.setAll(TB_START);
                    if (i == 1) {
                        Ik[i] = sc.gap_open;
                        tb.setIBits(TB_START);
                    } else {
                        const i_score = sc.gap_open + sc.gap_extend * (@as(i32, @intCast(i)) - 1);
                        const c_score = sc.xclip_prefix + sc.gap_open;
                        if (i_score > c_score) {
                            Ik[i] = i_score;
                            tb.setIBits(TB_INS);
                        } else {
                            Ik[i] = c_score;
                            tb.setIBits(TB_XCLIP_PREFIX);
                        }
                    }

                    if (i == m) {
                        tb.setSBits(TB_XCLIP_SUFFIX);
                    } else {
                        Sk[i] = MIN_SCORE;
                    }

                    if (Ik[i] > Sk[i]) {
                        Sk[i] = Ik[i];
                        tb.setSBits(TB_INS);
                    }
                    if (sc.xclip_prefix > Sk[i]) {
                        Sk[i] = sc.xclip_prefix;
                        tb.setSBits(TB_XCLIP_PREFIX);
                    }
                    if (i != m and Sk[i] + sc.xclip_suffix > Sk[m]) {
                        Sk[m] = Sk[i] + sc.xclip_suffix;
                        self.Lx.items[0] = m - i;
                    }
                    if (k == 0) self.traceback.set(i, 0, tb);
                    if (Sk[i] + sc.yclip_suffix > self.Sn.items[i]) {
                        self.Sn.items[i] = Sk[i] + sc.yclip_suffix;
                        self.Ly.items[i] = n;
                    }
                }
            }

            var j: usize = 1;
            while (j <= n) : (j += 1) {
                const curr = j % 2;
                const prev = 1 - curr;

                {
                    var tb = TracebackCell{};
                    self.I[curr].items[0] = MIN_SCORE;
                    if (j == 1) {
                        self.D[curr].items[0] = sc.gap_open;
                        tb.setDBits(TB_START);
                    } else {
                        const d_score = sc.gap_open + sc.gap_extend * (@as(i32, @intCast(j)) - 1);
                        const c_score = sc.yclip_prefix + sc.gap_open;
                        if (d_score > c_score) {
                            self.D[curr].items[0] = d_score;
                            tb.setDBits(TB_DEL);
                        } else {
                            self.D[curr].items[0] = c_score;
                            tb.setDBits(TB_YCLIP_PREFIX);
                        }
                    }
                    if (self.D[curr].items[0] > sc.yclip_prefix) {
                        self.S[curr].items[0] = self.D[curr].items[0];
                        tb.setSBits(TB_DEL);
                    } else {
                        self.S[curr].items[0] = sc.yclip_prefix;
                        tb.setSBits(TB_YCLIP_PREFIX);
                    }
                    if (j == n and self.Sn.items[0] > self.S[curr].items[0]) {
                        self.S[curr].items[0] = self.Sn.items[0];
                        tb.setSBits(TB_YCLIP_SUFFIX);
                    } else if (self.S[curr].items[0] + sc.yclip_suffix > self.Sn.items[0]) {
                        self.Sn.items[0] = self.S[curr].items[0] + sc.yclip_suffix;
                        self.Ly.items[0] = n - j;
                    }
                    self.traceback.set(0, j, tb);
                }

                {
                    var i: usize = 1;
                    while (i <= m) : (i += 1) self.S[curr].items[i] = MIN_SCORE;
                }

                const q = y[j - 1];
                const xclip_score = sc.xclip_prefix +
                    @max(sc.yclip_prefix, sc.gap_open + sc.gap_extend * (@as(i32, @intCast(j)) - 1));

                const Icurr = self.I[curr].items;
                const Dcurr = self.D[curr].items;
                const Dprev = self.D[prev].items;
                const Scurr = self.S[curr].items;
                const Sprev = self.S[prev].items;

                var i: usize = 1;
                while (i <= m) : (i += 1) {
                    const p = x[i - 1];
                    var tb = TracebackCell{};

                    const m_score = Sprev[i - 1] + sc.match_fn.score(p, q);

                    const i_score = Icurr[i - 1] + sc.gap_extend;
                    const s_score_i = Scurr[i - 1] + sc.gap_open;
                    var best_i_score: i32 = undefined;
                    if (i_score > s_score_i) {
                        best_i_score = i_score;
                        tb.setIBits(TB_INS);
                    } else {
                        best_i_score = s_score_i;
                        tb.setIBits(self.traceback.get(i - 1, j).getSBits());
                    }

                    const d_score = Dprev[i] + sc.gap_extend;
                    const s_score_d = Sprev[i] + sc.gap_open;
                    var best_d_score: i32 = undefined;
                    if (d_score > s_score_d) {
                        best_d_score = d_score;
                        tb.setDBits(TB_DEL);
                    } else {
                        best_d_score = s_score_d;
                        tb.setDBits(self.traceback.get(i, j - 1).getSBits());
                    }

                    tb.setSBits(TB_XCLIP_SUFFIX);
                    var best_s_score = Scurr[i];

                    if (m_score > best_s_score) {
                        best_s_score = m_score;
                        tb.setSBits(if (p == q) TB_MATCH else TB_SUBST);
                    }
                    if (best_i_score > best_s_score) {
                        best_s_score = best_i_score;
                        tb.setSBits(TB_INS);
                    }
                    if (best_d_score > best_s_score) {
                        best_s_score = best_d_score;
                        tb.setSBits(TB_DEL);
                    }
                    if (xclip_score > best_s_score) {
                        best_s_score = xclip_score;
                        tb.setSBits(TB_XCLIP_PREFIX);
                    }
                    const yclip_score = sc.yclip_prefix + sc.gap_open + sc.gap_extend * (@as(i32, @intCast(i)) - 1);
                    if (yclip_score > best_s_score) {
                        best_s_score = yclip_score;
                        tb.setSBits(TB_YCLIP_PREFIX);
                    }

                    Scurr[i] = best_s_score;
                    Icurr[i] = best_i_score;
                    Dcurr[i] = best_d_score;

                    if (Scurr[i] + sc.xclip_suffix > Scurr[m]) {
                        Scurr[m] = Scurr[i] + sc.xclip_suffix;
                        self.Lx.items[j] = m - i;
                    }
                    if (Scurr[i] + sc.yclip_suffix > self.Sn.items[i]) {
                        self.Sn.items[i] = Scurr[i] + sc.yclip_suffix;
                        self.Ly.items[i] = n - j;
                    }
                    self.traceback.set(i, j, tb);
                }
            }

            // Suffix clipping in the j = n column.
            {
                const jn = n;
                const curr = n % 2;
                const Scurr = self.S[curr].items;
                var i: usize = 0;
                while (i <= m) : (i += 1) {
                    if (self.Sn.items[i] > Scurr[i]) {
                        Scurr[i] = self.Sn.items[i];
                        self.traceback.getPtr(i, jn).setSBits(TB_YCLIP_SUFFIX);
                    }
                    if (Scurr[i] + sc.xclip_suffix > Scurr[m]) {
                        Scurr[m] = Scurr[i] + sc.xclip_suffix;
                        self.Lx.items[jn] = m - i;
                        self.traceback.getPtr(m, jn).setSBits(TB_XCLIP_SUFFIX);
                    }
                }
            }

            // Recompute the last column of I after S changed above.
            {
                const jn = n;
                const curr = n % 2;
                const Scurr = self.S[curr].items;
                const Icurr = self.I[curr].items;
                var i: usize = 1;
                while (i <= m) : (i += 1) {
                    const s_score = Scurr[i - 1] + sc.gap_open;
                    if (s_score > Icurr[i]) {
                        Icurr[i] = s_score;
                        const s_bit = self.traceback.get(i - 1, jn).getSBits();
                        self.traceback.getPtr(i, jn).setIBits(s_bit);
                    }
                    if (s_score > Scurr[i]) {
                        Scurr[i] = s_score;
                        self.traceback.getPtr(i, jn).setSBits(TB_INS);
                        if (Scurr[i] + sc.xclip_suffix > Scurr[m]) {
                            Scurr[m] = Scurr[i] + sc.xclip_suffix;
                            self.Lx.items[jn] = m - i;
                            self.traceback.getPtr(m, jn).setSBits(TB_XCLIP_SUFFIX);
                        }
                    }
                }
            }

            // Traceback.
            var i = m;
            var j2 = n;
            var operations: std.ArrayListUnmanaged(AlignmentOperation) = .empty;
            errdefer operations.deinit(self.gpa);
            var xstart: usize = 0;
            var ystart: usize = 0;
            var xend = m;
            var yend = n;

            var last_layer = self.traceback.get(i, j2).getSBits();
            while (true) {
                var next_layer: u16 = undefined;
                switch (last_layer) {
                    TB_START => break,
                    TB_INS => {
                        try operations.append(self.gpa, .ins);
                        next_layer = self.traceback.get(i, j2).getIBits();
                        i -= 1;
                    },
                    TB_DEL => {
                        try operations.append(self.gpa, .del);
                        next_layer = self.traceback.get(i, j2).getDBits();
                        j2 -= 1;
                    },
                    TB_MATCH => {
                        try operations.append(self.gpa, .match);
                        next_layer = self.traceback.get(i - 1, j2 - 1).getSBits();
                        i -= 1;
                        j2 -= 1;
                    },
                    TB_SUBST => {
                        try operations.append(self.gpa, .subst);
                        next_layer = self.traceback.get(i - 1, j2 - 1).getSBits();
                        i -= 1;
                        j2 -= 1;
                    },
                    TB_XCLIP_PREFIX => {
                        try operations.append(self.gpa, .{ .xclip = i });
                        xstart = i;
                        i = 0;
                        next_layer = self.traceback.get(0, j2).getSBits();
                    },
                    TB_XCLIP_SUFFIX => {
                        try operations.append(self.gpa, .{ .xclip = self.Lx.items[j2] });
                        i -= self.Lx.items[j2];
                        xend = i;
                        next_layer = self.traceback.get(i, j2).getSBits();
                    },
                    TB_YCLIP_PREFIX => {
                        try operations.append(self.gpa, .{ .yclip = j2 });
                        ystart = j2;
                        j2 = 0;
                        next_layer = self.traceback.get(i, 0).getSBits();
                    },
                    TB_YCLIP_SUFFIX => {
                        try operations.append(self.gpa, .{ .yclip = self.Ly.items[i] });
                        j2 -= self.Ly.items[i];
                        yend = j2;
                        next_layer = self.traceback.get(i, j2).getSBits();
                    },
                    else => unreachable,
                }
                last_layer = next_layer;
            }

            std.mem.reverse(AlignmentOperation, operations.items);
            const ops = try operations.toOwnedSlice(self.gpa);

            return .{
                .score = self.S[n % 2].items[m],
                .ystart = ystart,
                .xstart = xstart,
                .yend = yend,
                .xend = xend,
                .ylen = n,
                .xlen = m,
                .operations = ops,
                .mode = .custom,
            };
        }

        fn saveClips(self: *Self) [4]i32 {
            return .{ self.scoring.xclip_prefix, self.scoring.xclip_suffix, self.scoring.yclip_prefix, self.scoring.yclip_suffix };
        }
        fn restoreClips(self: *Self, c: [4]i32) void {
            self.scoring.xclip_prefix = c[0];
            self.scoring.xclip_suffix = c[1];
            self.scoring.yclip_prefix = c[2];
            self.scoring.yclip_suffix = c[3];
        }

        pub fn global(self: *Self, x: []const u8, y: []const u8) !Alignment {
            const saved = self.saveClips();
            self.scoring.xclip_prefix = MIN_SCORE;
            self.scoring.xclip_suffix = MIN_SCORE;
            self.scoring.yclip_prefix = MIN_SCORE;
            self.scoring.yclip_suffix = MIN_SCORE;
            var a = try self.custom(x, y);
            a.mode = .global;
            self.restoreClips(saved);
            return a;
        }

        pub fn semiglobal(self: *Self, x: []const u8, y: []const u8) !Alignment {
            const saved = self.saveClips();
            self.scoring.xclip_prefix = MIN_SCORE;
            self.scoring.xclip_suffix = MIN_SCORE;
            self.scoring.yclip_prefix = 0;
            self.scoring.yclip_suffix = 0;
            var a = try self.custom(x, y);
            a.mode = .semiglobal;
            try a.filterClipOperations(self.gpa);
            self.restoreClips(saved);
            return a;
        }

        pub fn local(self: *Self, x: []const u8, y: []const u8) !Alignment {
            const saved = self.saveClips();
            self.scoring.xclip_prefix = 0;
            self.scoring.xclip_suffix = 0;
            self.scoring.yclip_prefix = 0;
            self.scoring.yclip_suffix = 0;
            var a = try self.custom(x, y);
            a.mode = .local;
            try a.filterClipOperations(self.gpa);
            self.restoreClips(saved);
            return a;
        }
    };
}

// Tests -----------------------------------------------------------------------

const scores = @import("../scores.zig");

const Op = AlignmentOperation;

fn expectOps(expected: []const Op, actual: []const Op) !void {
    try testing.expectEqualDeep(expected, actual);
}

test "traceback cell packing" {
    var tb = TracebackCell{};
    tb.setAll(TB_SUBST);
    try testing.expectEqual(TB_SUBST, tb.getIBits());
    try testing.expectEqual(TB_SUBST, tb.getDBits());
    try testing.expectEqual(TB_SUBST, tb.getSBits());
    tb.setDBits(TB_INS);
    try testing.expectEqual(TB_INS, tb.getDBits());
    tb.setIBits(TB_XCLIP_PREFIX);
    try testing.expectEqual(TB_INS, tb.getDBits());
    try testing.expectEqual(TB_XCLIP_PREFIX, tb.getIBits());
    tb.setSBits(TB_YCLIP_SUFFIX);
    try testing.expectEqual(TB_YCLIP_SUFFIX, tb.getSBits());
    try testing.expectEqual(TB_XCLIP_PREFIX, tb.getIBits());
}

test "semiglobal" {
    var a = Aligner(MatchParams).init(testing.allocator, -5, -1, MatchParams.init(1, -1));
    defer a.deinit();
    var aln = try a.semiglobal("ACCGTGGAT", "AAAAACCGTTGAT");
    defer aln.deinit(testing.allocator);
    try testing.expectEqual(@as(usize, 4), aln.ystart);
    try testing.expectEqual(@as(usize, 0), aln.xstart);
    try expectOps(&.{ .match, .match, .match, .match, .match, .subst, .match, .match, .match }, aln.operations);
}

test "semiglobal gap_open < mismatch (score underflow guard)" {
    var a = Aligner(MatchParams).init(testing.allocator, -1, -1, MatchParams.init(1, -5));
    defer a.deinit();
    var aln = try a.semiglobal("ACCGTGGAT", "AAAAACCGTTGAT");
    defer aln.deinit(testing.allocator);
    try testing.expectEqual(@as(usize, 4), aln.ystart);
    try testing.expectEqual(@as(usize, 0), aln.xstart);
    try expectOps(&.{ .match, .match, .match, .match, .del, .match, .ins, .match, .match, .match }, aln.operations);
}

test "global" {
    var a = Aligner(MatchParams).init(testing.allocator, -5, -1, MatchParams.init(1, -1));
    defer a.deinit();
    var aln = try a.global("ACCGTGGAT", "AAAAACCGTTGAT");
    defer aln.deinit(testing.allocator);
    try testing.expectEqual(@as(usize, 0), aln.ystart);
    try testing.expectEqual(@as(usize, 0), aln.xstart);
    try expectOps(&.{ .del, .del, .del, .del, .match, .match, .match, .match, .match, .subst, .match, .match, .match }, aln.operations);
}

test "global affine ins" {
    var a = Aligner(MatchParams).init(testing.allocator, -5, -1, MatchParams.init(1, -3));
    defer a.deinit();
    var aln = try a.global("ACGAGAACA", "ACGACA");
    defer aln.deinit(testing.allocator);
    try expectOps(&.{ .match, .match, .match, .ins, .ins, .ins, .match, .match, .match }, aln.operations);
}

test "local" {
    var a = Aligner(MatchParams).init(testing.allocator, -5, -1, MatchParams.init(1, -1));
    defer a.deinit();
    var aln = try a.local("ACCGTGGAT", "AAAAACCGTTGAT");
    defer aln.deinit(testing.allocator);
    try testing.expectEqual(@as(usize, 4), aln.ystart);
    try testing.expectEqual(@as(usize, 0), aln.xstart);
    try expectOps(&.{ .match, .match, .match, .match, .match, .subst, .match, .match, .match }, aln.operations);
    try testing.expectEqual(@as(i32, 7), aln.score);
}

test "global blosum62" {
    var a = Aligner(FnMatch(scores.blosum62)).init(testing.allocator, -5, -1, .{});
    defer a.deinit();
    var aln = try a.global("AAAA", "AAAA");
    defer aln.deinit(testing.allocator);
    try testing.expectEqual(@as(i32, 16), aln.score);
    try expectOps(&.{ .match, .match, .match, .match }, aln.operations);
}

test "local blosum62" {
    var a = Aligner(FnMatch(scores.blosum62)).init(testing.allocator, -10, -1, .{});
    defer a.deinit();
    var aln = try a.local("LSPADKTNVKAA", "PEEKSAV");
    defer aln.deinit(testing.allocator);
    try testing.expectEqual(@as(usize, 2), aln.xstart);
    try testing.expectEqual(@as(usize, 9), aln.xend);
    try testing.expectEqual(@as(usize, 0), aln.ystart);
    try testing.expectEqual(@as(usize, 7), aln.yend);
    try expectOps(&.{ .match, .subst, .subst, .match, .subst, .subst, .match }, aln.operations);
    try testing.expectEqual(@as(i32, 16), aln.score);
}

test "custom semiglobal via clips (explicit Yclip)" {
    var sc = scoringFromScores(-5, -1, 1, -1).xclip(MIN_SCORE).yclip(0);
    var a = Aligner(MatchParams).withScoring(testing.allocator, sc);
    defer a.deinit();
    _ = &sc;
    var aln = try a.custom("ACCGTGGAT", "AAAAACCGTTGAT");
    defer aln.deinit(testing.allocator);
    try testing.expectEqual(@as(usize, 4), aln.ystart);
    try testing.expectEqual(@as(usize, 0), aln.xstart);
    try expectOps(&.{ .{ .yclip = 4 }, .match, .match, .match, .match, .match, .subst, .match, .match, .match }, aln.operations);
}
