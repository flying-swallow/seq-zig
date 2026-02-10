const std = @import("std");
const seq = @import("root.zig");

pub fn faWriteRecord(writer: *std.Io.Writer, slice: struct {
    name: []const u8,
    seq: []const u8,
}, options: struct {
    line_ending: enum {
        LF,
        CRLF,
    } = .LF,
}) !void {
    const line_ending = switch (options.line_ending) {
        .LF => "\n",
        .CRLF => "\r\n",
    };

    try writer.writeByte('>');
    try writer.writeAll(slice.name);
    try writer.writeAll(line_ending);

    try writer.writeAll(slice.seq);
    try writer.writeAll(line_ending);
}

pub fn fqWriteRecord(writer: *std.Io.Writer, slice: struct {
    name: []const u8,
    qual: []const u8,
    seq: []const u8,
}, options: struct {
    line_ending: enum {
        LF,
        CRLF,
    } = .LF,
}) !void {
    const line_ending = switch (options.line_ending) {
        .LF => "\n",
        .CRLF => "\r\n",
    };

    try writer.writeByte('@');
    try writer.writeAll(slice.name);
    try writer.writeAll(line_ending);

    try writer.writeAll(slice.seq);
    try writer.writeAll(line_ending);

    try writer.writeByte('+');
    try writer.writeAll(line_ending);

    try writer.writeAll(slice.qual);
    try writer.writeAll(line_ending);
}


pub fn takeFqSequence(reader: *std.Io.Reader, name: *std.Io.Writer, sequence: *std.Io.Writer, qual: *std.Io.Writer) !bool {
    const format: seq.Format = switch (try reader.takeByte()) {
        '@' => .fq,
        '>' => .fa,
        else => |c| {
            std.debug.print("Invalid first byte: {d}\n", .{c});
            return error.RecordParseError;
        },
    };
    try name.writeAll(std.mem.trimEnd(u8, try reader.takeDelimiterInclusive('\n'), &std.ascii.whitespace));

    var qual_line_count: usize = 0;

    var seq_line_count: usize = 0;
    var seq_number_bases: usize = 0;
    var qual_number_bases: usize = 0;

    var bases_per_line_expected: ?usize = null;

    var bases_per_line: usize = 0;

    while (true) {
        switch (reader.peekByte() catch |err| switch (err) {
            error.EndOfStream => return false, // end of stream no more records
            else => return err,
        }) {
            '+' => {
                _ = try reader.takeDelimiter('\n');
                while (true) {
                    if (qual_number_bases == seq_number_bases) {
                        if (seq_line_count != qual_line_count) return error.RecordParseError;
                        return true;
                    } else if (qual_number_bases > seq_number_bases) {
                        return error.RecordParseError;
                    }

                    if (qual_line_count > 0) {
                        if (bases_per_line_expected) |count| {
                            if (bases_per_line != count) return error.RecordParseError;
                        } else return error.RecordParseError;
                    }

                    qual_line_count += 1;

                    var qual_slice = std.mem.trimEnd(u8, try reader.takeDelimiterInclusive('\n'), &std.ascii.whitespace);
                    try qual.writeAll(qual_slice);

                    bases_per_line = qual_slice.len;
                    qual_number_bases += bases_per_line;
                }
            },
            '>' => {
                if (format == .fq)
                    return error.RecordParseError;
                // we finished reading this record
                return true;
            },
            else => {},
        }
        if (seq_line_count > 0) {
            if (bases_per_line_expected) |count| {
                if (bases_per_line != count) return error.RecordParseError;
            } else bases_per_line_expected = bases_per_line;
        }
        seq_line_count += 1;

        var seq_slice = std.mem.trimEnd(u8, try reader.takeDelimiterInclusive('\n'), &std.ascii.whitespace);
        try sequence.writeAll(seq_slice);

        bases_per_line = seq_slice.len;
        seq_number_bases += bases_per_line;
    }
}

test "read fq scanner" {
    const file = @embedFile("./test/t1.fq");
    const TestCase = struct { name: []const u8, sequence: []const u8, quality: []const u8 };
    var reader: std.Io.Reader = .fixed(file[0..]);
    //var fqReader: FqReader = .init(std.testing.allocator, &reader);
    //defer fqReader.deinit();

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

    var name: std.Io.Writer.Allocating = .init(std.testing.allocator);
    var sequence: std.Io.Writer.Allocating = .init(std.testing.allocator);
    var qual: std.Io.Writer.Allocating = .init(std.testing.allocator);
    defer name.deinit();
    defer sequence.deinit();
    defer qual.deinit();

    for (test_cases) |case| {
        try std.testing.expect(try takeFqSequence(&reader, &name.writer, &sequence.writer, &qual.writer));
        try std.testing.expectEqualStrings(case.name, name.written());
        try std.testing.expectEqualStrings(case.quality, qual.written());
        try std.testing.expectEqualStrings(case.sequence, sequence.written());
        
        name.clearRetainingCapacity();
        sequence.clearRetainingCapacity();
        qual.clearRetainingCapacity();
    }
    // zig fmt: on
}
