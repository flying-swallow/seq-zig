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
    pub const FQRecord = struct { name: []const u8 = undefined, seq: []const u8 = undefined, qual: []const u8 = undefined };

    arena: std.heap.ArenaAllocator,
    reader: *std.Io.Reader,

    pub fn init(allocator: std.mem.Allocator, reader: *std.Io.Reader) FQReader {
        return .{ .arena = .init(allocator), .reader = reader };
    }

    pub fn next(self: *FQReader) !FQRecord {
        self.arena.reset(.retain_capacity);
        const format: seq.SeqContianer = switch (try self.reader.takeByte()) {
            '@' => .fastq,
            '>' => .fasta,
            else => return error.RecordParseError,
        };
        var name_writer: std.Io.Writer.Allocating = .init(self.arena.allocator());
        var seq_writer: std.Io.Writer.Allocating = .init(self.arena.allocator());
        var qual_writer: std.Io.Writer.Allocating = .init(self.arena.allocator());
        try self.reader.streamDelimiter(&name_writer, '\n');

        var seq_line_count: u32 = 0;
        var qual_line_count: u32 = 0;

        var seq_number_bases: usize = 0;
        var qual_number_bases: usize = 0;

        //var bytes_per_line: usize = 0;
        var bytes_per_line_expected: ?usize = null;

        //var bases_per_line: usize = 0;
        var bases_per_line_expected: ?usize = null;

        var process_state: enum {
            seq,
            qual,
        } = .seq;

        while (true) finish: {
            switch (try self.reader.peekByte()) {
                '+' => {
                    _ = try self.reader.takeDelimiter('\n');
                    while (true) {
                        if (qual_number_bases == seq_number_bases) {
                            if (seq_line_count != qual_line_count) return error.RecordParseError;
                            break :finish;
                        } else if (self.qual_number_bases > self.seq_number_bases) {
                            return error.RecordParseError;
                        }
                        //if (qual_line_count > 0) {
                        //    if (bytes_per_line_expected) |count| {
                        //        if (bytes_per_line != count) return error.RecordParseError;
                        //    } else return error.RecordParseError;
                        //    if (bases_per_line_expected) |count| {
                        //        if (bases_per_line != count) return error.RecordParseError;
                        //    } else return error.RecordParseError;
                        //}
                    }
                },
                '>' => {
                    if (format == .fastq)
                        return error.RecordParseError;
                    break :finish;
                },
                _ => {},
            }
            seq_line_count += 1;
            const bytes_read = try self.reader.streamDelimiter(&seq_writer, '\n');
            var written_buf = seq_writer.written();
            const buffer_trimmed_len = std.mem.trimEnd(u8, written_buf, &std.ascii.whitespace).len;
            seq_writer.shrinkRetainingCapacity(buffer_trimmed_len);

            const bytes_per_line = bytes_read;
            const bases_per_line = bytes_read - (written_buf.len - buffer_trimmed_len);
            seq_number_bases += buffer_trimmed_len;

            if (seq_line_count > 0) {
                // we're onto the next line need to check line_counts
                if (bytes_per_line_expected) |count| {
                    if (bytes_per_line != count) return error.RecordParseError;
                } else bytes_per_line_expected = bytes_per_line;
                if (bases_per_line_expected) |count| {
                    if (bases_per_line != count) return error.RecordParseError;
                } else bases_per_line_expected = bases_per_line;
            }
        }

        self.current = .{
            .name = try name_writer.written(),
        };

        //while (true) {
        //    switch (self.state) {
        //        .start => {
        //            if (self.current) |rec| {
        //                self.current = null;
        //                return rec;
        //            } else {}
        //        },
        //        .name => {
        //            var name_writer: std.Io.Writer.Allocating = .init(self.arena.allocator());
        //            try self.reader.streamDelimiter(&name_writer, '\n');
        //            self.current = .{
        //                .name = try name_writer.toOwnedSlice(),
        //            };
        //            self.state = .seq;
        //        },
        //    }
        //}
        return .{};
    }
};
