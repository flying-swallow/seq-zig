# zbio

A small, dependency-free Zig library for reading and writing common bioinformatics
sequence formats.

> **Status:** early work in progress — the API is unstable and may change.

## Requirements

- Zig `0.17.0-dev` (master). zbio uses the new `std.Io.Reader` / `std.Io.Writer`
  interfaces, so it will not build on 0.13/0.14 releases.

## Modules

Everything is re-exported from the root module (`@import("zbio")`):

| Module      | What it does |
| ----------- | ------------ |
| `fa`        | FASTA / FASTQ record reading and writing |
| `fai`       | `.fai` index parsing |
| `sam`       | SAM text reader — header + alignment records (CIGAR, optional fields, FLAG) |
| `bam`       | BAM reader **and** writer, including the BGZF block codec |
| `alphabet`  | Sequence alphabets (DNA/RNA/protein), complement / reverse-complement, IUPAC `Nucleotide` & `AminoAcid` enums, FASTQ quality encodings, and `RankTransform` / q-grams |
| `stats`     | Length statistics (N50/L50, quartiles) and GC content / counts |
| `bitenc`    | Fixed-width (1–8 bit) sequence bit-packing |
| `scores`    | BLOSUM (30/45/62) and PAM (40/120/200/250) substitution matrices |
| `orf`       | Forward-strand open reading frame finder |
| `alignment` | `distance` (Hamming / Levenshtein / bounded) and `pairwise` (generalized Smith-Waterman: global / semiglobal / local, affine gaps) |
| `pm`        | Exact pattern matching — Horspool, KMP, Shift-And, BNDM, BOM |
| `index`     | Suffix array (+ LCP), BWT (`less` / `Occ`), and FM-index backward search |

The `alignment`, `pm`, and `index` modules were ported from
[rust-bio](https://github.com/rust-bio/rust-bio). Suffix arrays are built by
prefix doubling (rather than rust-bio's SAIS); the result is identical since the
suffix array is unique. Deferred for now: Myers approximate matching, the
bidirectional FMD-index (SMEMs), and the sampled suffix array.

## Use it in your project

Fetch it (point at wherever this repo lives):

```sh
zig fetch --save git+https://github.com/you/zbio.git
```

Wire it up in `build.zig`:

```zig
const zbio = b.dependency("zbio", .{ .target = target, .optimize = optimize });
exe.root_module.addImport("zbio", zbio.module("zbio"));
```

and import it:

```zig
const zbio = @import("zbio");
```

## Examples

### Read FASTA / FASTQ records

```zig
const std = @import("std");
const zbio = @import("zbio");

var name: std.Io.Writer.Allocating = .init(allocator);
var seq: std.Io.Writer.Allocating = .init(allocator);
var qual: std.Io.Writer.Allocating = .init(allocator);
defer { name.deinit(); seq.deinit(); qual.deinit(); }

var r: std.Io.Reader = .fixed(fastq_bytes);
while (try zbio.fa.takeFqSequence(&r, &name.writer, &seq.writer, &qual.writer)) {
    std.debug.print(">{s}\n{s}\n", .{ name.written(), seq.written() });
    name.clearRetainingCapacity();
    seq.clearRetainingCapacity();
    qual.clearRetainingCapacity();
}
```

### Read a SAM file

```zig
const sam = zbio.sam;

var r: std.Io.Reader = .fixed(sam_bytes);
var sr = sam.reader(allocator, &r);
var header = try sr.readHeader();
defer header.deinit(allocator);

while (try sr.next(&header)) |rec| {
    std.debug.print("{s} {s}:{d} flag={d}\n", .{
        rec.qname, header.refName(rec.tid) orelse "*", rec.pos, rec.flag,
    });
    var cigar = rec.cigar();
    while (try cigar.next()) |op| {
        std.debug.print("  {d}{c}\n", .{ op.len, op.op.toChar() });
    }
    var aux = rec.aux();
    while (try aux.next()) |tag| { _ = tag; }
}
```

### Read a BAM file (BGZF is decompressed for you)

```zig
const bam = zbio.bam;

var r: std.Io.Reader = .fixed(bam_bytes);
var reader = try bam.open(allocator, &r);
defer reader.deinit();

while (try reader.next()) |rec| {
    const flags = rec.flags();               // typed FLAG bits
    if (flags.unmap) continue;
    const name = reader.header.refName(rec.ref_id) orelse "*";
    std.debug.print("{s} {s}:{d}\n", .{ rec.qname, name, rec.pos });
    // rec.seqBase(i), rec.cigar(), rec.aux() decode the packed fields
}
```

### Write a BAM file

```zig
var out: std.Io.Writer.Allocating = .init(allocator);
defer out.deinit();

var writer = bam.writer(allocator, &out.writer);
defer writer.deinit();
try writer.writeHeader(reader.header_text, reader.header.refs);
while (try reader.next()) |rec| try writer.writeRecord(&rec);
try writer.finish();                          // flushes the final block + BGZF EOF marker
// out.written() now holds a valid .bam byte stream
```

### Alphabets, complement, and quality

```zig
const ab = zbio.alphabet;

_ = ab.dna.isValid('A');           // true (matching is case-insensitive)
_ = ab.dna.complement('a');        // 't'

var buf: [4]u8 = undefined;
ab.dna.reverseComplementInto(&buf, "AACG");        // buf == "CGTT"

_ = ab.Nucleotide.fromChar('R');                   // IUPAC ambiguity code (4-bit)
_ = ab.AminoAcid.fromChar('W').?.threeLetter();    // "Trp"
_ = ab.Encoding.sanger.decodeToPhred('I');         // 40
```

## Build & test

```sh
zig build          # build the library + test artifact
zig build test     # run the unit tests
zig build bench    # run the benchmark suite (ReleaseFast)
```

## Benchmarks

Benchmarks use [zBench](https://github.com/hendriknielaender/zBench). It is a
**dev-only, lazy dependency** — declared with `.lazy = true` in `build.zig.zon`,
so it is fetched only when you run `zig build bench`. The published `zbio`
library module stays dependency-free; `zig build test` never compiles or links
zBench. The bench executable is always built in `ReleaseFast`, independent of
`-Doptimize`, and reports throughput (MB/s, items/s) for both the format-IO code
and the algorithm modules. Benchmark sources live in `bench/`.

## Not yet implemented

- BAI / CSI indexing and region queries
- CRAM
- VCF / BCF, BED
- Streaming (non-in-memory) BGZF decode for very large BAMs

## References

- SAM/BAM specification: <https://samtools.github.io/hts-specs/SAMv1.pdf>
- The `alphabet` module is modeled on
  [biogo/alphabet](https://github.com/biogo/biogo/tree/master/alphabet).
