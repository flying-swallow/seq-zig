const std = @import("std");
const seq = @import("root.zig");

const testing = std.testing;

// http://www.htslib.org/doc/faidx.html
// faidx – an index enabling random access to FASTA and FASTQ files
// if the reserved is 0 then we assume its a const []u8 else the entry takes up some amount of reserved memory
pub const FaiFaIndex = struct {
    name: ?[]const u8 = null,
    length: u32 = 0, // Total length of this reference sequence, in bases
    offset: u64 = 0, // Offset in the FASTA/FASTQ file of this sequence's first base
    line_base: u32 = 0, // The number of bases on each line
    line_width: u32 = 0, // The number of bytes in each line, including the newline
};

pub const FaiFqIndex = struct {
    name: ?[]const u8 = null,
    length: u32 = 0, // Total length of this reference sequence, in bases
    offset: u64 = 0, // Offset in the FASTA/FASTQ file of this sequence's first base
    line_base: u32 = 0, // The number of bases on each line
    line_width: u32 = 0, // The number of bytes in each line, including the newline
    qual_offset: u64 = 0, // Offset of sequence's first quality within the FASTQ file

    pub fn read(self: *FaiFqIndex, reader: *std.Io.File.Reader, sequence: *std.Io.Writer, qual: *std.Io.Writer) !void {
        try reader.seekTo(self.offset);
        var remaining_bases: usize = self.length;
        while (remaining_bases > 0) {
            const bases_to_read: usize = @min(self.line_base, remaining_bases);
            sequence.writeAll(try reader.interface.take(bases_to_read));
            remaining_bases -= bases_to_read;
            if (remaining_bases > 0)
                try reader.interface.toss(self.line_width - self.line_base);
        }

        try reader.seekTo(self.qual_offset);
        remaining_bases = self.length;
        while (remaining_bases > 0) {
            const bases_to_read: usize = @min(self.line_base, remaining_bases);
            qual.writeAll(try reader.interface.take(bases_to_read));
            remaining_bases -= bases_to_read;
            if (remaining_bases > 0)
                try reader.interface.toss(self.line_width - self.line_base);
        }
    }
};

// The sequence name recorded in a .fai index is the record header up to the
// first whitespace; anything after it (a description) is ignored for indexing.
fn faiNameToken(header: []const u8) []const u8 {
    const trimmed = std.mem.trimEnd(u8, header, &std.ascii.whitespace);
    return if (std.mem.indexOfAny(u8, trimmed, &std.ascii.whitespace)) |i| trimmed[0..i] else trimmed;
}

// scan a reader and produce a fasta index for random access
// into a fasta file
pub const FaiFaIndexIndexScanner = struct {
    offset: usize = 0,
    reader: *std.Io.Reader,
    allocator: std.mem.Allocator,
    name_buf: std.ArrayList(u8) = .empty,

    pub fn init(allocator: std.mem.Allocator, reader: *std.Io.Reader) FaiFaIndexIndexScanner {
        return .{ .allocator = allocator, .reader = reader };
    }

    pub fn deinit(self: *FaiFaIndexIndexScanner) void {
        self.name_buf.deinit(self.allocator);
    }

    // Scan the next FASTA record and produce its .fai entry, or null at end of
    // stream. The returned `name` borrows the scanner's buffer until the next
    // call.
    pub fn next(self: *FaiFaIndexIndexScanner) !?FaiFaIndex {
        const marker = self.reader.takeByte() catch |err| switch (err) {
            error.EndOfStream => return null,
            else => |e| return e,
        };
        if (marker != '>') return error.RecordParseError;
        self.offset += 1;

        const header = try self.reader.takeDelimiterInclusive('\n');
        self.offset += header.len;
        self.name_buf.clearRetainingCapacity();
        try self.name_buf.appendSlice(self.allocator, faiNameToken(header));

        const seq_offset = self.offset;

        var length: u32 = 0;
        var line_base: u32 = 0;
        var line_width: u32 = 0;
        var line_count: usize = 0;
        var bases_per_line: usize = 0;
        var bases_per_line_expected: ?usize = null;

        while (true) {
            switch (self.reader.peekByte() catch |err| switch (err) {
                error.EndOfStream => break, // end of stream, last record read
                else => |e| return e,
            }) {
                '>' => break, // start of the next record
                else => {},
            }

            // The previous line now has a successor, so every line but the last
            // of a record must share the same width.
            if (line_count > 0) {
                if (bases_per_line_expected) |count| {
                    if (bases_per_line != count) return error.RecordParseError;
                } else bases_per_line_expected = bases_per_line;
            }
            line_count += 1;

            const line = try self.reader.takeDelimiterInclusive('\n');
            self.offset += line.len;
            bases_per_line = std.mem.trimEnd(u8, line, &std.ascii.whitespace).len;

            if (line_count == 1) {
                line_base = @intCast(bases_per_line);
                line_width = @intCast(line.len);
            }
            length += @intCast(bases_per_line);
        }

        return .{
            .name = self.name_buf.items,
            .length = length,
            .offset = seq_offset,
            .line_base = line_base,
            .line_width = line_width,
        };
    }
};

// scan a reader and produce a fastq index for random access
// into a fastq file
pub const FaiFqIndexIndexScanner = struct {
    offset: usize = 0,
    reader: *std.Io.Reader,
    allocator: std.mem.Allocator,
    name_buf: std.ArrayList(u8) = .empty,

    pub fn init(allocator: std.mem.Allocator, reader: *std.Io.Reader) FaiFqIndexIndexScanner {
        return .{ .allocator = allocator, .reader = reader };
    }

    pub fn deinit(self: *FaiFqIndexIndexScanner) void {
        self.name_buf.deinit(self.allocator);
    }

    // Scan the next FASTQ record and produce its .fai entry, or null at end of
    // stream. The returned `name` borrows the scanner's buffer until the next
    // call. Quality lines are consumed by base count (not by marker) because a
    // quality line may itself begin with '@' or '+'.
    pub fn next(self: *FaiFqIndexIndexScanner) !?FaiFqIndex {
        const marker = self.reader.takeByte() catch |err| switch (err) {
            error.EndOfStream => return null,
            else => |e| return e,
        };
        if (marker != '@') return error.RecordParseError;
        self.offset += 1;

        const header = try self.reader.takeDelimiterInclusive('\n');
        self.offset += header.len;
        self.name_buf.clearRetainingCapacity();
        try self.name_buf.appendSlice(self.allocator, faiNameToken(header));

        const seq_offset = self.offset;

        var length: u32 = 0;
        var line_base: u32 = 0;
        var line_width: u32 = 0;

        var seq_line_count: usize = 0;
        var qual_line_count: usize = 0;
        var seq_number_bases: usize = 0;
        var qual_number_bases: usize = 0;
        var bases_per_line: usize = 0;
        var bases_per_line_expected: ?usize = null;

        while (true) {
            switch (self.reader.peekByte() catch |err| switch (err) {
                error.EndOfStream => return error.RecordParseError, // truncated record: no quality
                else => |e| return e,
            }) {
                '+' => {
                    // Consume the '+' separator line; quality starts after it.
                    const sep = try self.reader.takeDelimiterInclusive('\n');
                    self.offset += sep.len;
                    const qual_offset = self.offset;

                    while (true) {
                        if (qual_number_bases == seq_number_bases) {
                            if (seq_line_count != qual_line_count) return error.RecordParseError;
                            return .{
                                .name = self.name_buf.items,
                                .length = length,
                                .offset = seq_offset,
                                .line_base = line_base,
                                .line_width = line_width,
                                .qual_offset = qual_offset,
                            };
                        } else if (qual_number_bases > seq_number_bases) {
                            return error.RecordParseError;
                        }

                        if (qual_line_count > 0) {
                            if (bases_per_line_expected) |count| {
                                if (bases_per_line != count) return error.RecordParseError;
                            } else return error.RecordParseError;
                        }
                        qual_line_count += 1;

                        const qual_line = try self.reader.takeDelimiterInclusive('\n');
                        self.offset += qual_line.len;
                        bases_per_line = std.mem.trimEnd(u8, qual_line, &std.ascii.whitespace).len;
                        qual_number_bases += bases_per_line;
                    }
                },
                else => {},
            }

            // Validate the previous sequence line's width now that a successor
            // exists; the final sequence line may be shorter.
            if (seq_line_count > 0) {
                if (bases_per_line_expected) |count| {
                    if (bases_per_line != count) return error.RecordParseError;
                } else bases_per_line_expected = bases_per_line;
            }
            seq_line_count += 1;

            const line = try self.reader.takeDelimiterInclusive('\n');
            self.offset += line.len;
            bases_per_line = std.mem.trimEnd(u8, line, &std.ascii.whitespace).len;

            if (seq_line_count == 1) {
                line_base = @intCast(bases_per_line);
                line_width = @intCast(line.len);
            }
            seq_number_bases += bases_per_line;
            length += @intCast(bases_per_line);
        }
    }
};

pub fn takeFaFaiIndex(reader: *std.Io.Reader) !?FaiFaIndex {
    var col: u8 = 0;
    var length: u32 = 0;
    var offset: u64 = 0;
    var line_base: u32 = 0;
    var line_width: u32 = 0;
    var name: []const u8 = &.{};

    while (true) {
        var buf = reader.takeDelimiterInclusive('\n') catch |err| switch (err) {
            error.EndOfStream => return null,
            else => |e| return e,
        };
        while (buf.len > 0) {
            const column = std.mem.trimEnd(u8, if (std.mem.findScalar(u8, buf, std.ascii.control_code.ht)) |end| res: {
                const result = buf[0..end];
                buf = buf[(end + 1)..];
                break :res result;
            } else res: {
                const result = buf[0..];
                buf = &.{};
                break :res result;
            }, &std.ascii.whitespace);
            switch (col) {
                0 => name = std.mem.trimEnd(u8, column, &std.ascii.whitespace),
                1 => length = std.fmt.parseUnsigned(u32, column, 10) catch return error.ParseError,
                2 => offset = std.fmt.parseUnsigned(u64, column, 10) catch return error.ParseError,
                3 => line_base = std.fmt.parseUnsigned(u32, column, 10) catch return error.ParseError,
                4 => line_width = std.fmt.parseUnsigned(u32, column, 10) catch return error.ParseError,
                else => return error.UnexpectedNumberColum,
            }
            col += 1;
        }
        if (col > 0) {
            if (col != 5)
                return error.UnexpectedNumberColum; // the column
            return .{
                .name = name,
                .length = length,
                .offset = offset,
                .line_base = line_base,
                .line_width = line_width,
            };
        }
    }
}

pub fn takeFqFaiIndex(reader: *std.Io.Reader) !?FaiFqIndex {
    var col: u8 = 0;
    var length: u32 = 0;
    var offset: u64 = 0;
    var line_base: u32 = 0;
    var line_width: u32 = 0;
    var qual_offset: u64 = 0;
    var name: []const u8 = &.{};

    while (true) {
        var buf = reader.takeDelimiterInclusive('\n') catch |err| switch (err) {
            error.EndOfStream => return null,
            else => |e| return e,
        };
        while (buf.len > 0) {
            const column = std.mem.trimEnd(u8, if (std.mem.findScalar(u8, buf, std.ascii.control_code.ht)) |end| res: {
                const result = buf[0..end];
                buf = buf[(end + 1)..];
                break :res result;
            } else res: {
                const result = buf[0..];
                buf = &.{};
                break :res result;
            }, &std.ascii.whitespace);
            switch (col) {
                0 => name = std.mem.trimEnd(u8, column, &std.ascii.whitespace),
                1 => length = std.fmt.parseUnsigned(u32, column, 10) catch return error.ParseError,
                2 => offset = std.fmt.parseUnsigned(u64, column, 10) catch return error.ParseError,
                3 => line_base = std.fmt.parseUnsigned(u32, column, 10) catch return error.ParseError,
                4 => line_width = std.fmt.parseUnsigned(u32, column, 10) catch return error.ParseError,
                5 => qual_offset = std.fmt.parseUnsigned(u64, column, 10) catch return error.ParseError,
                else => return error.UnexpectedNumberColum,
            }
            col += 1;
        }
        if (col > 0) {
            if (col != 6)
                return error.UnexpectedNumberColum; // the column
            return .{ .name = name, .length = length, .offset = offset, .line_base = line_base, .line_width = line_width, .qual_offset = qual_offset };
        }
    }
}

pub const FaqIndexReader = struct {
    reader: *std.Io.Reader,

    pub fn init(reader: *std.Io.Reader) FaqIndexReader {
        return .{ .reader = reader };
    }

    pub fn next(self: *FaqIndexReader) !?FaiFqIndex {
        var col: u8 = 0;
        var length: u32 = 0;
        var offset: u64 = 0;
        var line_base: u32 = 0;
        var line_width: u32 = 0;
        var qual_offset: u64 = 0;
        var name: []const u8 = &.{};

        while (true) {
            var buf = self.reader.takeDelimiterInclusive('\n') catch |err| switch (err) {
                error.EndOfStream => return null,
                else => |e| return e,
            };
            while (buf.len > 0) {
                const column = std.mem.trimEnd(u8, if (std.mem.findScalar(u8, buf, std.ascii.control_code.ht)) |end| res: {
                    const result = buf[0..end];
                    buf = buf[(end + 1)..];
                    break :res result;
                } else res: {
                    const result = buf[0..];
                    buf = &.{};
                    break :res result;
                }, &std.ascii.whitespace);
                switch (col) {
                    0 => name = std.mem.trimEnd(u8, column, &std.ascii.whitespace),
                    1 => length = std.fmt.parseUnsigned(u32, column, 10) catch return error.ParseError,
                    2 => offset = std.fmt.parseUnsigned(u64, column, 10) catch return error.ParseError,
                    3 => line_base = std.fmt.parseUnsigned(u32, column, 10) catch return error.ParseError,
                    4 => line_width = std.fmt.parseUnsigned(u32, column, 10) catch return error.ParseError,
                    5 => qual_offset = std.fmt.parseUnsigned(u64, column, 10) catch return error.ParseError,
                    else => return error.UnexpectedNumberColum,
                }
                col += 1;
            }
            if (col > 0) {
                if (col != 6)
                    return error.UnexpectedNumberColum; // the column
                return .{ .name = name, .length = length, .offset = offset, .line_base = line_base, .line_width = line_width, .qual_offset = qual_offset };
            }
        }
    }
};

test "fasta parse fai index" {
    const file = @embedFile("./test/ce.fa.fai");
    var reader: std.Io.Reader = .fixed(file[0..]);

    const test_cases =
        [_]struct {
            name: []const u8,
            length: u32,
            offset: u32,
            line_base: u32,
            line_width: u32,
        }{ .{ .name = "CHROMOSOME_I", .length = 1009800, .offset = 14, .line_base = 50, .line_width = 51 }, .{ .name = "CHROMOSOME_II", .length = 5000, .offset = 1030025, .line_base = 50, .line_width = 51 }, .{ .name = "CHROMOSOME_III", .length = 5000, .offset = 1035141, .line_base = 50, .line_width = 51 }, .{ .name = "CHROMOSOME_IV", .length = 5000, .offset = 1040256, .line_base = 50, .line_width = 51 }, .{ .name = "CHROMOSOME_V", .length = 5000, .offset = 1045370, .line_base = 50, .line_width = 51 }, .{ .name = "CHROMOSOME_X", .length = 5000, .offset = 1050484, .line_base = 50, .line_width = 51 }, .{ .name = "CHROMOSOME_MtDNA", .length = 5000, .offset = 1055602, .line_base = 50, .line_width = 51 } };

    for (test_cases) |case| {
        if (try takeFaFaiIndex(&reader)) |record| {
            try testing.expectEqualStrings(case.name, record.name.?);
            try testing.expectEqual(@as(u64, case.length), record.length);
            try testing.expectEqual(@as(u64, case.offset), record.offset);
            try testing.expectEqual(@as(u32, case.line_base), record.line_base);
            try testing.expectEqual(@as(u32, case.line_width), record.line_width);
        } else try testing.expect(false);
    }
}

test "fastq parse fq.fai index" {
    const file = @embedFile("./test/t3.fq.fai");
    var reader: std.Io.Reader = .fixed(file[0..]);
    //var iter: FaqIndexReader = .init(&reader);
    const test_cases =
        [_]struct { name: []const u8, length: u32, offset: u32, line_base: u32, line_width: u32, qual_offset: u32 }{
            .{ .name = "FAKE0005_1", .length = 63, .offset = 85, .line_base = 63, .line_width = 64, .qual_offset = 151 },
            .{ .name = "FAKE0006_1", .length = 63, .offset = 300, .line_base = 63, .line_width = 64, .qual_offset = 366 },
            .{ .name = "FAKE0005_2", .length = 63, .offset = 515, .line_base = 63, .line_width = 64, .qual_offset = 581 },
            .{ .name = "FAKE0006_2", .length = 63, .offset = 730, .line_base = 63, .line_width = 64, .qual_offset = 796 },
            .{ .name = "FAKE0005_3", .length = 63, .offset = 945, .line_base = 63, .line_width = 64, .qual_offset = 1011 },
            .{ .name = "FAKE0006_3", .length = 63, .offset = 1160, .line_base = 63, .line_width = 64, .qual_offset = 1226 },
            .{ .name = "FAKE0005_4", .length = 63, .offset = 1375, .line_base = 63, .line_width = 64, .qual_offset = 1441 },
            .{ .name = "FAKE0006_4", .length = 63, .offset = 1590, .line_base = 63, .line_width = 64, .qual_offset = 1656 },
        };
    for (test_cases) |case| {
        if (try takeFqFaiIndex(&reader)) |entry| {
            try testing.expectEqualStrings(case.name, entry.name.?);
            try testing.expectEqual(@as(u64, case.length), entry.length);
            try testing.expectEqual(@as(u64, case.offset), entry.offset);
            try testing.expectEqual(@as(u32, case.line_base), entry.line_base);
            try testing.expectEqual(@as(u32, case.line_width), entry.line_width);
            try testing.expectEqual(@as(u32, case.qual_offset), entry.qual_offset);
        } else try testing.expect(false);
    }
}

test "fasta scan fai index" {
    // Scan the FASTA data file and verify the scanner reproduces the same index
    // that `samtools faidx` produced (src/test/ce.fa.fai).
    const file = @embedFile("./test/ce.fa");
    var reader: std.Io.Reader = .fixed(file[0..]);
    var scanner: FaiFaIndexIndexScanner = .init(testing.allocator, &reader);
    defer scanner.deinit();

    const test_cases =
        [_]struct {
            name: []const u8,
            length: u32,
            offset: u32,
            line_base: u32,
            line_width: u32,
        }{ .{ .name = "CHROMOSOME_I", .length = 1009800, .offset = 14, .line_base = 50, .line_width = 51 }, .{ .name = "CHROMOSOME_II", .length = 5000, .offset = 1030025, .line_base = 50, .line_width = 51 }, .{ .name = "CHROMOSOME_III", .length = 5000, .offset = 1035141, .line_base = 50, .line_width = 51 }, .{ .name = "CHROMOSOME_IV", .length = 5000, .offset = 1040256, .line_base = 50, .line_width = 51 }, .{ .name = "CHROMOSOME_V", .length = 5000, .offset = 1045370, .line_base = 50, .line_width = 51 }, .{ .name = "CHROMOSOME_X", .length = 5000, .offset = 1050484, .line_base = 50, .line_width = 51 }, .{ .name = "CHROMOSOME_MtDNA", .length = 5000, .offset = 1055602, .line_base = 50, .line_width = 51 } };

    for (test_cases) |case| {
        if (try scanner.next()) |record| {
            try testing.expectEqualStrings(case.name, record.name.?);
            try testing.expectEqual(@as(u64, case.length), record.length);
            try testing.expectEqual(@as(u64, case.offset), record.offset);
            try testing.expectEqual(@as(u32, case.line_base), record.line_base);
            try testing.expectEqual(@as(u32, case.line_width), record.line_width);
        } else try testing.expect(false);
    }
    try testing.expect((try scanner.next()) == null);
}

test "fastq scan fq.fai index" {
    // Scan the FASTQ data file and verify the scanner reproduces the same index
    // that `samtools faidx` produced (src/test/t3.fq.fai). t3.fq holds 105
    // records; we check the leading ones, which also exercise a quality line
    // that begins with '@' (t3.fq line 4).
    const file = @embedFile("./test/t3.fq");
    var reader: std.Io.Reader = .fixed(file[0..]);
    var scanner: FaiFqIndexIndexScanner = .init(testing.allocator, &reader);
    defer scanner.deinit();

    const test_cases =
        [_]struct { name: []const u8, length: u32, offset: u32, line_base: u32, line_width: u32, qual_offset: u32 }{
            .{ .name = "FAKE0005_1", .length = 63, .offset = 85, .line_base = 63, .line_width = 64, .qual_offset = 151 },
            .{ .name = "FAKE0006_1", .length = 63, .offset = 300, .line_base = 63, .line_width = 64, .qual_offset = 366 },
            .{ .name = "FAKE0005_2", .length = 63, .offset = 515, .line_base = 63, .line_width = 64, .qual_offset = 581 },
            .{ .name = "FAKE0006_2", .length = 63, .offset = 730, .line_base = 63, .line_width = 64, .qual_offset = 796 },
            .{ .name = "FAKE0005_3", .length = 63, .offset = 945, .line_base = 63, .line_width = 64, .qual_offset = 1011 },
            .{ .name = "FAKE0006_3", .length = 63, .offset = 1160, .line_base = 63, .line_width = 64, .qual_offset = 1226 },
            .{ .name = "FAKE0005_4", .length = 63, .offset = 1375, .line_base = 63, .line_width = 64, .qual_offset = 1441 },
            .{ .name = "FAKE0006_4", .length = 63, .offset = 1590, .line_base = 63, .line_width = 64, .qual_offset = 1656 },
        };
    for (test_cases) |case| {
        if (try scanner.next()) |entry| {
            try testing.expectEqualStrings(case.name, entry.name.?);
            try testing.expectEqual(@as(u64, case.length), entry.length);
            try testing.expectEqual(@as(u64, case.offset), entry.offset);
            try testing.expectEqual(@as(u32, case.line_base), entry.line_base);
            try testing.expectEqual(@as(u32, case.line_width), entry.line_width);
            try testing.expectEqual(@as(u32, case.qual_offset), entry.qual_offset);
        } else try testing.expect(false);
    }
}
