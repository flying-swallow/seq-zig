//! Benchmark/test corpora embedded at build time. Rooted in src/test/ so these
//! @embedFile paths stay inside this module's own root directory (embedding them
//! from bench/ would escape the module root and fail to compile). Imported only
//! by the bench executable, so `zig build test` never compiles this data.

pub const ce_fa: []const u8 = @embedFile("ce.fa"); // ~1 MB multi-line FASTA
pub const ce_fa_fai: []const u8 = @embedFile("ce.fa.fai"); // its .fai index
pub const t3_fq: []const u8 = @embedFile("t3.fq"); // multi-record FASTQ
pub const t3_fq_fai: []const u8 = @embedFile("t3.fq.fai");
pub const t1_fq: []const u8 = @embedFile("t1.fq");
pub const t1_bam: []const u8 = @embedFile("t1.bam"); // BGZF BAM
