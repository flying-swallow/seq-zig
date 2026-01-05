const seq = @import("./seq.zig");
const std = @import("std");
const testing = std.testing;

pub const FqWriteOption = struct {
    lineEnding: enum { crlf, lf } = .lf,
};

pub fn FqWriter(comptime Stream: type) type {
    return struct {
        const Self = @This();

        stream: *Stream,
        options: FqWriteOption,
        pub const WriterError = Stream.Writer.Error;

        pub fn write_record(self: *Self, slice: struct {
            name: []const u8,
            qual: []const u8,
            seq: []const u8,
        }) Stream.Writer.Error!void {
            var fq = self.stream.writer();

            try fq.writeByte('@');
            try fq.writeAll(slice.name);
            try fq.writeByte('\n');

            try fq.writeAll(slice.seq);
            try fq.writeByte('\n');

            try fq.writeByte('+');
            try fq.writeByte('\n');

            try fq.writeAll(slice.qual);
            try fq.writeByte('\n');
        }
    };
}

pub fn fqWriter(stream: anytype, option: FqWriteOption) FqWriter(@TypeOf(stream.*)) {
    return .{
        .options = option,
        .stream = stream,
    };
}

pub const ReadIterOptions = struct {
    reserve_len: usize = 2048,
};

fn streamIntoArrayUntilDelim(stream: anytype, buffer: []u8, delimter: u8) error{StreamTooLong, EndOfStream}![]u8{
    var pos: usize= 0;
    while(true) {
        const byte = try stream.readByte();
        if(byte == delimter) break;
        buffer[pos] = byte;
        pos += 1;
        if(pos == buffer.len) return error.StreamTooLong;
    }
    return buffer[0..pos];
}



pub fn Scanner(comptime Stream: type, comptime options: ReadIterOptions) type {
    return struct {
        const Self = @This();

        stream: Stream,
        buffer: [options.reserve_len]u8 = undefined, // staged data on the stack for the iter
        pos: usize = 0,
        state: enum { start, name, seq, seq_end, qual, qual_end} = .start,
        format: seq.SeqContainer = .fastq,

        seq_line_count: u32 = 0,
        qual_line_count: u32 = 0,

        seq_number_bases: usize = 0,
        qual_number_bases: usize = 0,

        bytes_per_line: usize = 0,
        bytes_per_line_expected: ?usize = null,

        bases_per_line: usize = 0,
        bases_per_line_expected: ?usize = null,

        buffer_exausted: bool = false,

        fn internalReset(self: *Self) void {
            self.seq_line_count = 0;
            self.qual_line_count = 0;
            self.seq_number_bases = 0;
            self.qual_number_bases = 0;
            self.bytes_per_line = 0;
            self.bytes_per_line_expected = null;
            self.bases_per_line = 0;
            self.bases_per_line_expected = null;
            self.buffer_exausted = false;
        }

        pub fn next(self: *Self) (Stream.ReadError || error{ BufferOverflow, RecordParseError })!?union(enum) {
            begin: struct {buf: []const u8},
            seq: struct {finished: bool, buf: []const u8},
            qual: struct {finished: bool, buf: []const u8},
            end: struct {
             //   length: u32 = 0, // Total length of this reference sequence, in bases
             //   offset: u64 = 0, // Offset in the FASTA/FASTQ file of this sequence's first base
             //   line_base: u32 = 0, // The number of bases on each line
             //   line_width: u32 = 0, // The number of bytes in each line, including the newline
             //   qual_offset: ?u64 = null, // Offset of sequence's first quality within the FASTQ file
            },
        } {
            var reader = self.stream.reader();

            self.pos = 0;
            while (true) {
                switch (self.state) {
                    .start => {
                        const byte: u8 = reader.readByte() catch |err| switch (err) {
                            error.EndOfStream => return null, // we just have a name so the buffer is malformed
                            else => |e| return e,
                        };
                        self.state = switch (byte) {
                            '@' => l: {
                                self.format = .fastq;
                                break :l .name;
                            },
                            '>' => l: {
                                self.format = .fasta;
                                break :l .name;
                            },
                            else => return error.RecordParseError,
                        };
                    },
                    .name => {
                        var bufferStream = std.io.fixedBufferStream(&self.buffer);
                        reader.streamUntilDelimiter(bufferStream.writer(), '\n', null) catch |err| switch (err) {
                            error.EndOfStream => return null,
                            error.NoSpaceLeft, error.StreamTooLong => return error.BufferOverflow,
                            else => |e| return e,
                        };
                        self.pos = @as(usize, bufferStream.pos);
                        self.state = .seq;
                        return .{ .begin = .{
                            .buf = self.buffer[0..self.pos]
                        }};
                    },
                    .seq_end => {
                        self.internalReset();
                        self.state = .name;
                        return .{ .end = .{}};
                    },
                    .seq =>  {
                        if (!self.buffer_exausted) {
                            switch (reader.readByte() catch |err| switch (err) {
                                error.EndOfStream => return null,
                            }) {
                                '+' => {
                                    self.bases_per_line = 0;
                                    self.bytes_per_line = 0;
                                    self.state = .qual;
                                    try reader.skipUntilDelimiterOrEof('\n');
                                    return .{ .seq = .{
                                        .finished = true,
                                        .buf = self.buffer[0..self.pos]
                                    }};

                                },
                                '>' => {
                                    if (self.format == .fastq)
                                        return error.RecordParseError;
                                    self.internalReset();
                                    self.state = .seq_end;
                                    return .{ .seq = .{.finished = true, .buf = self.buffer[0..self.pos]}};
                                },
                                else => |b| {
                                    if (self.seq_line_count > 0) {
                                        // we're onto the next line need to check line_counts
                                        if (self.bytes_per_line_expected) |count| {
                                            if (self.bytes_per_line != count) return error.RecordParseError;
                                        } else self.bytes_per_line_expected = self.bytes_per_line;
                                        if (self.bases_per_line_expected) |count| {
                                            if (self.bases_per_line != count) return error.RecordParseError;
                                        } else self.bases_per_line_expected = self.bases_per_line;
                                    }
                                    self.bases_per_line = 0;
                                    self.bytes_per_line = 0;
                                    self.buffer[self.pos] = b;
                                    self.pos += 1;
                                    self.bytes_per_line += 1;
                                    self.bases_per_line += 1;
                                    self.seq_number_bases += 1;
                                },
                            }
                        }

                        self.seq_line_count += 1;
                        self.buffer_exausted = false;

                        var byte_buf = streamIntoArrayUntilDelim(reader, self.buffer[self.pos..], '\n') catch |err| switch(err) {
                            error.EndOfStream => return null,
                            error.StreamTooLong => blk: {
                                self.seq_line_count -= 1;
                                self.buffer_exausted = true;
                                break :blk self.buffer[self.pos..];
                            },
                            else => |e| return e,
                        };
                        const bases_buf = std.mem.trimRight(u8, byte_buf, &std.ascii.whitespace);

                        self.bytes_per_line += byte_buf.len;
                        self.bases_per_line += bases_buf.len;
                        self.seq_number_bases += bases_buf.len;
                        self.pos += bases_buf.len;

                        if (self.buffer_exausted) {
                            return .{ .seq = .{
                                .finished = false,
                                .buf = self.buffer[0..self.pos]
                            }};
                        }
                    },
                    .qual_end => {
                        self.internalReset();
                        self.state = .start;
                        return .{ .end = .{}};
                    },
                    .qual => {
                        if (!self.buffer_exausted) {
                            if (self.qual_number_bases == self.seq_number_bases) {
                                if (self.seq_line_count != self.qual_line_count) return error.RecordParseError;
                                self.state = .qual_end;
                                return .{ .qual = .{.finished = true, .buf = self.buffer[0..self.pos]}};
                            } else if (self.qual_number_bases > self.seq_number_bases) {
                                return error.RecordParseError;
                            }
                            if (self.qual_line_count > 0) {
                                if (self.bytes_per_line_expected) |count| {
                                    if (self.bytes_per_line != count) return error.RecordParseError;
                                } else return error.RecordParseError;
                                if (self.bases_per_line_expected) |count| {
                                    if (self.bases_per_line != count) return error.RecordParseError;
                                } else return error.RecordParseError;
                            }
                            self.bases_per_line = 0;
                            self.bytes_per_line = 0;
                        }

                        self.qual_line_count += 1;
                        self.buffer_exausted = false;
                        
                        var byte_buf = streamIntoArrayUntilDelim(reader, self.buffer[self.pos..], '\n') catch |err| switch(err) {
                            error.EndOfStream => return null,
                            error.StreamTooLong => blk: {
                                self.seq_line_count -= 1;
                                self.buffer_exausted = true;
                                break :blk self.buffer[self.pos..];
                            },
                            else => |e| return e,
                        };
                        const bases_buf = std.mem.trimRight(u8, byte_buf, &std.ascii.whitespace);

                        self.bytes_per_line += byte_buf.len;
                        self.bases_per_line += bases_buf.len;
                        self.qual_number_bases += bases_buf.len;
                        self.pos += bases_buf.len;
                        
                        if (self.buffer_exausted) {
                            return .{ .qual = .{.finished = false, .buf = self.buffer[0..self.pos]}};
                        }
                    },
                }
            }
            unreachable;
        }
    };
}

pub fn scanner(stream: anytype, comptime options: ReadIterOptions) Scanner(@TypeOf(stream), options) {
    return .{ .stream = stream };
}


test "read fq scanner" {
    const file = @embedFile("./test/t1.fq");
    var stream = std.io.fixedBufferStream(file);
    var scan = seq.fa.scanner(stream, .{});
    const TestCase = struct { name: []const u8, sequence: []const u8, quality: []const u8 };

    // zig fmt: off
    const test_cases =
        [_]TestCase{ 
        .{ 
            .name = "HWI-D00523:240:HF3WGBCXX:1:1101:2574:2226 1:N:0:CTGTAG", 
            .sequence = "TGAGGAATATTGGTCAATGGGCGCGAGCCTGAACCAGCCAAGTAGCGTGAAGGATGACTGCCCTACGGGTTGTAAACTTCTTTTATAAAGGAATAAAGTGAGGCACGTGTGCCTTTTTGTATGTACTTTATGAATAAGGATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATCCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGCAGGCGGT", 
            .quality = "HIHIIIIIHIIHGHHIHHIIIIIIIIIIIIIIIHHIIIIIHHIHIIIIIGIHIIIIHHHHHHGHIHIIIIIIIIIIIGHIIIIIGHIIIIHIIHIHHIIIIHIHHIIIIIIIGIIIIIIIHIIIIIGHIIIIHIIIH?DGHEEGHIIIIIIIIIIIHIIHIIIHHIIHIHHIHCHHIIHGIHHHHHHH<GG?B@EHDE-BEHHHII5B@GHHF?CGEHHHDHIHIIH" }, 
        .{ 
            .name = "HWI-D00523:240:HF3WGBCXX:1:1101:5586:3020 1:N:0:CTGTAG", 
            .sequence = "TGGGGAATATTGGGCAATGGGCGGAAGCCTGACCCAGCAACGCCGCGTGAAGGAAGAAGGCCCTCGGGTTGTAAACTTCTTTTCTATAGGACGAAGAAGTGACGGTACTATAGGAATAAGCCACGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGTGGCGAGCGTTATCCGGATTTACTGGGTGTAAAGGGCGTGTAGGCGGGAGAGCAAGTCAGATGTGA", 
            .quality = "EHEHGIIIHIGGHGHFEHHEHGCHHGGHIIIGHHIFHHGHHIEHIIIIGHHHHIIGIHGHGGHHHHHCHHHICCHHHHH@HHIIIGCEHGHHGHCHHGDGCGCCEHEEHGIIGHHGHHHIGGFCFHHIHIGIIIHGGHFHIIFEFHIIHIGDHCHFHHGCHCE?GHIIH<C?GHHHHIGFDEHHHHEC88@<@@EHHHIHDH-@HHCDHHDDEHH6@F6@6@@EH@@" }, 
        .{ 
            .name = "HWI-D00523:240:HF3WGBCXX:1:1101:2860:2149 1:N:0:CTGTAG", 
            .sequence = "TAGGGAATATTGCTCAATGGGGGAAACCCTGAAGCAGCAACGCCGCGTGGAGGATGAAGGTTTTAGGATTGTAAACTCCTTTTGTGAGAGAAGATTATGACGGTATCTCACGAATAAGCTCCGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGGAGCGAGCGTTGTCCGGAATTACTGGGTGTAAAGGGAGCGTAGGCGGGACTGCAAGTTGGGTGTCAAA", 
            .quality = "HGHHGHIIHIIIIIIHHHIIHIIIIIHGHCHIHIIHIIHIIIIIIIIIHIGHIEHIHIIG<FEHHHIHHIIIIIIHIFFHHHHIIIIIHHHHGHFEHHIHHHIHEHHFHHHIHIIIIIIIHIHDHHHHHEHHIIGIIHHIIGHHHIIGDDAGGHHFHHHIHICHHHGH,GHHHHGCEHEG?@6@G?-@>HHHHHHHDEH<@H-@CDD>:E?@GHEF-@E:@H+@-@@" 
        } 
    };
    // zig fmt: on

    var name = std.ArrayList(u8).init(std.testing.allocator);
    var sequence = std.ArrayList(u8).init(std.testing.allocator);
    var quality = std.ArrayList(u8).init(std.testing.allocator);
    defer {
        name.deinit();
        sequence.deinit();
        quality.deinit();
    }

    for (test_cases) |case| {
        blk: while (try scan.next()) |sec| {
            switch (sec) {
                .begin => |value| try name.appendSlice(value.buf),
                .seq => |value| try sequence.appendSlice(value.buf),
                .qual => |value| try quality.appendSlice(value.buf),
                .end => break :blk,
            }
        }
        try testing.expectEqualStrings(case.name, name.items);
        try testing.expectEqualStrings(case.quality, quality.items);
        try testing.expectEqualStrings(case.sequence, sequence.items);
  
        name.clearRetainingCapacity();
        sequence.clearRetainingCapacity();
        quality.clearRetainingCapacity();
    }
}
