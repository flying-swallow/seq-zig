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

    pub fn read(self: *FaiFqIndex, reader: *std.Io.File.Reader, sequence: *std.Io.Writer, qual: *std.Io.Writer) void {
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

// scan a reader and produce a fasta index for random access
// into a fasta file
pub const FaiFaIndexIndexScanner = struct {
    offset: usize = 0,
    reader: *std.Io.Reader,

    pub fn init(reader: *std.Io.Reader) FaiFaIndexIndexScanner {
        return .{ .reader = reader };
    }
};

pub const FaiFqIndexIndexScanner = struct {
    offset: usize = 0,
    reader: *std.Io.Reader,

    pub fn init(reader: *std.Io.Reader) FaiFaIndexIndexScanner {
        return .{ .reader = reader };
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
