const std = @import("std");
const seq = @import("root.zig");

pub const FQWriterOptions = struct {
    line_ending: enum {
        LF,
        CRLF,
    } = .LF,
};

pub const FQWriter = struct {
    writer: *std.Io.Writer,
    options: FQWriterOptions = .{},
    pub fn init(writer: *std.Io.Writer, options: FQWriterOptions) FQWriter {
        return .{ .writer = writer, .options = options };
    }

    pub fn writeRecord(self: *FQWriter, slice: struct {
        name: []const u8,
        qual: []const u8,
        seq: []const u8,
    }) !void {
        const line_ending = switch (self.options.line_ending) {
            .LF => "\n",
            .CRLF => "\r\n",
        };

        try self.writer.writeByte('@');
        try self.writer.writeAll(slice.name);
        try self.writer.writeAll(line_ending);

        try self.writer.writeAll(slice.seq);
        try self.writer.writeAll(line_ending);

        try self.writer.writeByte('+');
        try self.writer.writeAll(line_ending);

        try self.writer.writeAll(slice.qual);
        try self.writer.writeAll(line_ending);
    }
};

pub const FQReader = struct {
    pub const FQRecord = struct { name: []const u8 = &.{}, seq: []const u8 = &.{}, qual: []const u8 = &.{} };
    arena: std.heap.ArenaAllocator,
    reader: *std.Io.Reader,

    pub fn init(allocator: std.mem.Allocator, reader: *std.Io.Reader) FQReader {
        return .{ .arena = .init(allocator), .reader = reader };
    }

    pub fn deinit(self: *FQReader) void {
        self.arena.deinit();
    }

    pub fn next(self: *FQReader) !?FQRecord {
        _ = self.arena.reset(.retain_capacity);
        const format: seq.SeqContianer = switch (try self.reader.takeByte()) {
            '@' => .fastq,
            '>' => .fasta,
            else => |c| {
                std.debug.print("Invalid first byte: {d}\n", .{c});
                return error.RecordParseError;
            },
        };
        var name_writer: std.Io.Writer.Allocating = .init(self.arena.allocator());
        var seq_writer: std.Io.Writer.Allocating = .init(self.arena.allocator());
        var qual_writer: std.Io.Writer.Allocating = .init(self.arena.allocator());
        _ = try self.reader.streamDelimiterEnding(&name_writer.writer, '\n');
        self.reader.toss(1);

        var qual_line_count: usize = 0;

        var seq_line_count: usize = 0;
        var seq_number_bases: usize = 0;
        var qual_number_bases: usize = 0;

        var bases_per_line_expected: ?usize = null;

        var bases_per_line: usize = 0;

        while (true) {
            switch (self.reader.peekByte() catch |err| switch (err) {
                error.EndOfStream => return null,
                else => return err,
            }) {
                '+' => {
                    _ = try self.reader.takeDelimiter('\n');
                    while (true) {
                        if (qual_number_bases == seq_number_bases) {
                            if (seq_line_count != qual_line_count) return error.RecordParseError;
                            return .{
                                .name = name_writer.written(),
                                .seq = seq_writer.written(),
                                .qual = qual_writer.written(),
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

                        const last_len = qual_writer.written().len;
                        _ = try self.reader.streamDelimiterEnding(&qual_writer.writer, '\n');
                        self.reader.toss(1);
                        qual_writer.shrinkRetainingCapacity(std.mem.trimEnd(u8, qual_writer.written(), &std.ascii.whitespace).len);

                        bases_per_line = (qual_writer.written().len - last_len);
                        qual_number_bases += bases_per_line;
                    }
                },
                '>' => {
                    if (format == .fastq)
                        return error.RecordParseError;
                    return .{
                        .name = name_writer.written(),
                        .seq = seq_writer.written(),
                        .qual = qual_writer.written(),
                    };
                },
                else => {},
            }
            if (seq_line_count > 0) {
                if (bases_per_line_expected) |count| {
                    if (bases_per_line != count) return error.RecordParseError;
                } else bases_per_line_expected = bases_per_line;
            }
            seq_line_count += 1;

            const prev_written_len = seq_writer.written().len;
            _ = try self.reader.streamDelimiterEnding(&seq_writer.writer, '\n');
            self.reader.toss(1);
            seq_writer.shrinkRetainingCapacity(std.mem.trimEnd(u8, seq_writer.written(), &std.ascii.whitespace).len);

            bases_per_line = (seq_writer.written().len - prev_written_len);
            seq_number_bases += bases_per_line;
        }
    }
};

test "read fq scanner" {
    const file = @embedFile("./test/t1.fq");
    const TestCase = struct { name: []const u8, sequence: []const u8, quality: []const u8 };
    var reader: std.Io.Reader = .fixed(file[0..]);
    var fqReader: FQReader = .init(std.testing.allocator, &reader);
    defer fqReader.deinit();

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

    for (test_cases) |case| {
        var record = try fqReader.next();
        try std.testing.expectEqualStrings(case.name, record.?.name);
        try std.testing.expectEqualStrings(case.quality, record.?.qual);
        try std.testing.expectEqualStrings(case.sequence, record.?.seq);
    }
    // zig fmt: on
}
