const fa = @import("fa.zig");


pub const SeqContianer = enum {
    unknown,
    fastq, //  text-based format for nucleotide sequences 
    fasta,
    bam,

    fn name(con: SeqContianer) [] const u8 {
        switch(con) {
            .fastq => "FASTQ",
            .fasta => "FASTA",
            .bam => "BAM"
        }
    }

    pub const fastq_ext= [_]u8{"fq", "fasta"};
    pub const fasta_ext = [_]u8{"fa", "fastq"};

    fn ext(con: SeqContianer) [][] const u8 {
        switch (con) {
            .fasta => fasta_ext,
            .fastq => fastq_ext,
        }
        return {};
    }
};



pub const FQWriter = fa.FQWriter;
pub const FQReader = fa.FQReader;

test {
    
}
