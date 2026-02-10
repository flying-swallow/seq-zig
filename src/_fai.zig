const std = @import("std");
const seq = @import("seq.zig");
const testing = std.testing;

// zig fmt: off
const ColumnIdx = enum(u8) { 
    name = 0, 
    length = 1, 
    offset = 2, 
    line_bases = 3,
    line_width = 4, 
    qual_offset = 5 
};
// zig fmt: on

// http://www.htslib.org/doc/faidx.html
// faidx – an index enabling random access to FASTA and FASTQ files
// if the reserved is 0 then we assume its a const []u8 else the entry takes up some amount of reserved memory
pub const FaiFaIndex = struct {
    name: []const u8,
    length: u32 = 0, // Total length of this reference sequence, in bases
    offset: u64 = 0, // Offset in the FASTA/FASTQ file of this sequence's first base
    line_base: u32 = 0, // The number of bases on each line
    line_width: u32 = 0, // The number of bytes in each line, including the newline
};

pub const FaiFaqIndex = struct {
    name: []const u8,
    length: u32 = 0, // Total length of this reference sequence, in bases
    offset: u64 = 0, // Offset in the FASTA/FASTQ file of this sequence's first base
    line_base: u32 = 0, // The number of bases on each line
    line_width: u32 = 0, // The number of bytes in each line, including the newline
    qual_offset: u64 = 0, // Offset of sequence's first quality within the FASTQ file
};

pub const FaiReadIterOptions = struct {
    reserve_len: usize = 4096,
};

pub fn FastaIndexReadIterator(comptime Stream: type, comptime container: seq.SeqContainer, comptime options: FaiReadIterOptions) type {
    return struct {
        stream: Stream,
        buffer: [options.reserve_len]u8 = undefined,
        const Self = @This();

        pub fn next(self: *Self) (anyerror || error{ UnexpectedNumberColum, Overflow, ParseError })!?switch (container) {
            .fasta => FaiFaIndex,
            .fastq => FaiFaqIndex,
            else => unreachable,
        } {
            var reader = self.stream.reader();
            var col: u8 = 0;
            var pos: usize = 0;

            var name_len: usize = 0; // reserved the first chunk of the buffer for the name
            var length: u32 = 0;
            var offset: u64 = 0;
            var line_base: u32 = 0;
            var line_width: u32 = 0;
            var qual_offset: u64 = 0;

            while (true) {
                const byte: u8 = reader.readByte() catch |err| switch (err) {
                    error.EndOfStream => return null, // the buffer is exausted
                    else => |e| return e,
                };
                switch (byte) {
                    // skip \r
                    std.ascii.control_code.cr => break,
                    // horizontal tab
                    std.ascii.control_code.ht, std.ascii.control_code.lf => l: {
                        // line feed \n
                        if (std.ascii.control_code.lf == byte and col == 0) break :l; // if number of col is 0 then this is an empty row skip
                        switch (comptime container) {
                            .fastq => {
                                switch (@as(ColumnIdx, @enumFromInt(col))) {
                                    ColumnIdx.name => name_len = pos,
                                    ColumnIdx.length => length = std.fmt.parseUnsigned(u32, self.buffer[name_len..pos], 10) catch return error.ParseError,
                                    ColumnIdx.offset => offset = std.fmt.parseUnsigned(u64, self.buffer[name_len..pos], 10) catch return error.ParseError,
                                    ColumnIdx.line_bases => line_base = std.fmt.parseUnsigned(u32, self.buffer[name_len..pos], 10) catch return error.ParseError,
                                    ColumnIdx.line_width => line_width = std.fmt.parseUnsigned(u32, self.buffer[name_len..pos], 10) catch return error.ParseError,
                                    ColumnIdx.qual_offset => qual_offset = std.fmt.parseUnsigned(u64, self.buffer[name_len..pos], 10) catch return error.ParseError,
                                }
                            },
                            .fasta => {
                                switch (@as(ColumnIdx, @enumFromInt(col))) {
                                    ColumnIdx.name => name_len = pos,
                                    ColumnIdx.length => length = std.fmt.parseUnsigned(u32, self.buffer[name_len..pos], 10) catch return error.ParseError,
                                    ColumnIdx.offset => offset = std.fmt.parseUnsigned(u64, self.buffer[name_len..pos], 10) catch return error.ParseError,
                                    ColumnIdx.line_bases => line_base = std.fmt.parseUnsigned(u32, self.buffer[name_len..pos], 10) catch return error.ParseError,
                                    ColumnIdx.line_width => line_width = std.fmt.parseUnsigned(u32, self.buffer[name_len..pos], 10) catch return error.ParseError,
                                    else => return error.UnexpectedNumberColum,
                                }
                            },
                            else => unreachable,
                        }
                        pos = name_len;
                        col += 1;
                        if (std.ascii.control_code.lf == byte) break; // we terminate

                    },
                    else => {
                        self.buffer[pos] = byte;
                        pos += 1;
                    },
                }
            }
            switch (comptime container) {
                .fastq => {
                    if (col != 6)
                        return error.UnexpectedNumberColum; // the column
                    return .{ .name = self.buffer[0..name_len], .length = length, .offset = offset, .line_base = line_base, .line_width = line_width, .qual_offset = qual_offset };
                },
                .fasta => {
                    if (col != 5)
                        return error.UnexpectedNumberColum; // the column
                    return .{
                        .name = self.buffer[0..name_len],
                        .length = length,
                        .offset = offset,
                        .line_base = line_base,
                        .line_width = line_width,
                    };
                },
                else => unreachable,
            }
        }
    };
}

pub fn faiReaderIterator(stream: anytype, comptime container: seq.SeqContainer, comptime options: FaiReadIterOptions) FastaIndexReadIterator(@TypeOf(stream), container, options) {
    return .{ .stream = stream };
}

pub const FastaIndexOptions = struct {
    reserve_len: usize = 1024,
};

pub fn FastaIndexReader(comptime Stream: type, comptime container: seq.SeqContainer, comptime options: FastaIndexOptions) type {
    return struct {
        const Self = @This();

        stream: Stream,
        buffer: [options.reserve_len]u8 = undefined,

        pub fn read(
            self: *Self,
            index: switch (container) {
                .fasta => struct {
                    length: u32 = 0, // Total length of this reference sequence, in bases
                    offset: u64 = 0, // Offset in the FASTA/FASTQ file of this sequence's first base
                    line_base: u32 = 0, // The number of bases on each line
                    line_width: u32 = 0, // The number of bytes in each line, including the newline
                },
                .fastq => struct {
                    length: u32 = 0, // Total length of this reference sequence, in bases
                    offset: u64 = 0, // Offset in the FASTA/FASTQ file of this sequence's first base
                    line_base: u32 = 0, // The number of bases on each line
                    line_width: u32 = 0, // The number of bytes in each line, including the newline
                    qual_offset: u64 = 0, // Offset of sequence's first quality within the FASTQ file
                },
                else => unreachable,
            },
            record: anytype,
        ) (@TypeOf(record.*).Error || Stream.ReadError || error{ EndOfStream, BufferOverflow, RecordParseError })!bool {
            var reader = self.stream.reader();
            var seek = self.stream.seekableStream();
            {
                try seek.seekTo(index.offset); // seek to the offset of the fai
                var remaining_bases: usize = index.length;
                while (remaining_bases > 0) {
                    var bases_to_read: usize = @min(index.line_base, remaining_bases);
                    more_bases: while (true) {
                        const num_bytes: usize = try reader.read(self.buffer[0..@min(bases_to_read, options.reserve_len)]);
                        if (num_bytes == 0) return false;
                        try record.appendSeq(self.buffer[0..num_bytes]);
                        remaining_bases -= num_bytes;
                        bases_to_read -= num_bytes;

                        if (bases_to_read > 0) continue :more_bases;
                        break :more_bases;
                    }
                    try reader.skipBytes(index.line_width - index.line_base, .{ .buf_size = 16 });
                }
            }
            switch (container) {
                .fastq => {
                    try seek.seekTo(index.qual_offset); // seek to the offset of the fai
                    var remaining_bases: usize = index.length;
                    while (remaining_bases > 0) {
                        var bases_to_read: usize = @min(index.line_base, remaining_bases);
                        more_bases: while (true) {
                            const num_bytes: usize = try reader.read(self.buffer[0..@min(bases_to_read, options.reserve_len)]);
                            if (num_bytes == 0) return false;
                            try record.appendQual(self.buffer[0..num_bytes]);
                            remaining_bases -= num_bytes;
                            bases_to_read -= num_bytes;
                            if (bases_to_read > 0) continue :more_bases;
                            break :more_bases;
                        }
                        try reader.skipBytes(index.line_width - index.line_base, .{ .buf_size = 16 });
                    }
                },
                else => {},
            }
            return true;
        }
    };
}

pub fn fastaIndexReader(stream: anytype, comptime container: seq.SeqContainer, comptime options: FastaIndexOptions) FastaIndexReader(@TypeOf(stream), container, options) {
    return .{ .stream = stream };
}

test "parse fai index and seek fq" {
    //const indexFile = @embedFile("./test/t3.fq.fai");
    const fq = @embedFile("./test/t3.fq");
    var fqStream = std.io.fixedBufferStream(fq);
    const TestCase = struct { length: u32, offset: u64, line_base: u32, line_width: u32, qual_offset: u64, seq: []const u8, qual: []const u8 };

    const FaTestRecord = struct {
        const Self = @This();

        qual: std.ArrayListUnmanaged(u8),
        seq: std.ArrayListUnmanaged(u8),
        pub const Error = error{OutOfMemory};

        pub inline fn appendSeq(self: *Self, buf: []const u8) error{OutOfMemory}!void {
            try self.seq.appendSlice(std.testing.allocator, buf);
        }

        pub inline fn appendQual(self: *Self, buf: []const u8) error{OutOfMemory}!void {
            try self.qual.appendSlice(std.testing.allocator, buf);
        }
    };

    var faTestRecord: FaTestRecord = .{ .seq = .{}, .qual = .{} };
    defer {
        faTestRecord.qual.deinit(std.testing.allocator);
        faTestRecord.seq.deinit(std.testing.allocator);
    }

    var reader = fastaIndexReader(fqStream, .fastq, .{});
    // zig fmt: off
    const test_cases =
        [_]TestCase{
        .{ 
            .length = 63,  .offset = 85, .line_base = 63, .line_width = 64, .qual_offset = 151, 
            .seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACG", 
            .qual = "@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~" 
        },
        .{ 
            .length = 63, .offset = 1590, .line_base = 63, .line_width = 64, .qual_offset = 1656, 
            .seq = "GCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCA", 
            .qual = "~}|{zyxwvutsrqponmlkjihgfedcba`_^]\\[ZYXWVUTSRQPONMLKJIHGFEDCBA@" 
        },
    };
    // zig fmt: on

    for (test_cases) |case| {
        if (try reader.read(.{ .length = case.length, .offset = case.offset, .line_base = case.line_base, .line_width = case.line_width, .qual_offset = case.qual_offset }, &faTestRecord)) {
            try testing.expectEqualStrings(case.seq, faTestRecord.seq.items);
            try testing.expectEqualStrings(case.qual, faTestRecord.qual.items);
        } else try testing.expect(false);
        faTestRecord.qual.clearRetainingCapacity();
        faTestRecord.seq.clearRetainingCapacity();
    }
}

test "fastq parse fai index" {
    const file = @embedFile("./test/t3.fq.fai");
    var stream = std.io.fixedBufferStream(file);
    var reader = faiReaderIterator(stream, .fastq, .{});
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
        if (try reader.next()) |entry| {
            try testing.expectEqualStrings(case.name, entry.name);
            try testing.expectEqual(@as(u64, case.length), entry.length);
            try testing.expectEqual(@as(u64, case.offset), entry.offset);
            try testing.expectEqual(@as(u32, case.line_base), entry.line_base);
            try testing.expectEqual(@as(u32, case.line_width), entry.line_width);
        } else try testing.expect(false);
    }
}

test "fasta parse fai index" {
    const file = @embedFile("./test/ce.fa.fai");
    var stream = std.io.fixedBufferStream(file);
    var reader = faiReaderIterator(stream, .fasta, .{});

    const test_cases =
        [_]struct {
        name: []const u8,
        length: u32,
        offset: u32,
        line_base: u32,
        line_width: u32,
    }{ .{ .name = "CHROMOSOME_I", .length = 1009800, .offset = 14, .line_base = 50, .line_width = 51 }, .{ .name = "CHROMOSOME_II", .length = 5000, .offset = 1030025, .line_base = 50, .line_width = 51 }, .{ .name = "CHROMOSOME_III", .length = 5000, .offset = 1035141, .line_base = 50, .line_width = 51 }, .{ .name = "CHROMOSOME_IV", .length = 5000, .offset = 1040256, .line_base = 50, .line_width = 51 }, .{ .name = "CHROMOSOME_V", .length = 5000, .offset = 1045370, .line_base = 50, .line_width = 51 }, .{ .name = "CHROMOSOME_X", .length = 5000, .offset = 1050484, .line_base = 50, .line_width = 51 }, .{ .name = "CHROMOSOME_MtDNA", .length = 5000, .offset = 1055602, .line_base = 50, .line_width = 51 } };

    for (test_cases) |case| {
        if (try reader.next()) |record| {
            try testing.expectEqualStrings(case.name, record.name);
            try testing.expectEqual(@as(u64, case.length), record.length);
            try testing.expectEqual(@as(u64, case.offset), record.offset);
            try testing.expectEqual(@as(u32, case.line_base), record.line_base);
            try testing.expectEqual(@as(u32, case.line_width), record.line_width);
        } else try testing.expect(false);
    }
}
