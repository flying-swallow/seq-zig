pub const Format = enum {
    unknown,
    fa,
    fq, //  text-based format for nucleotide sequences
    fai,
    fqi,
    bam,
    bai,
    cram,
    crai,
    vcf,
    bcf,
    csi,
    gzi,
    tbi,
    bed,

    fn name(con: Format) []const u8 {
        switch (con) {
            .fq => "FASTQ",
            .fa => "FASTA",
            .bam => "BAM",
        }
    }

    pub const fastq_ext = [_]u8{ "fq", "fasta" };
    pub const fasta_ext = [_]u8{ "fa", "fastq" };

    fn ext(con: Format) [][]const u8 {
        switch (con) {
            .fa => fasta_ext,
            .fq => fastq_ext,
        }
        return {};
    }
};

pub const fa = @import("fa.zig");
pub const fai = @import("fai.zig");
pub const sam = @import("sam.zig");
pub const bam = @import("bam.zig");
pub const alphabet = @import("alphabet.zig");

test {
    _ = fa;
    _ = fai;
    _ = sam;
    _ = bam;
    _ = alphabet;
}
