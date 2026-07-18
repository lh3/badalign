# badalign

Extract *bad* read alignments from a BAM file. Two kinds of problems are reported,
both as TAB-delimited lines:

- **Type-1 (`C`)** — long clips (default ≥40 bp) at the primary alignment's clip
  boundaries. Useful for finding chimeric/misassembled reads and breakpoints.
- **Type-2 (`D`)** — mismatch-dense regions (default ≥5 high-quality ≥Q20
  mismatches within a 100 bp window). Useful for finding noisy / mis-mapped stretches.

Built on [`noodles`](https://github.com/zaeleus/noodles) for BAM/FASTA I/O.

## Build

```sh
cargo build --release
# binary at target/release/badalign
```

## Usage

```
badalign [options] aln.bam

Options:
  -r FILE   reference genome FASTA (loaded fully into RAM)
  -l INT    length of flanking regions around a clip / dense region [250]
  -Q INT    min base quality for a mismatch to count (type-2) [20]
  -w INT    window length to scan for a dense region (type-2) [100]
  -c INT    min clip length to report (type-1) [40]
  -m INT    min mismatches in a window (type-2) [5]
```

Output goes to stdout.

The input need not be sorted in any particular way — records are processed
streaming, and a read's chimeric segments are reconstructed from the primary
record's `SA` tag, so coordinate-sorted BAMs work as well as query-grouped ones.

## Output format

All coordinates are **0-based, half-open** (BED convention). Read-based positions
(`readStart/End`, `extractStart/End`, `clipOffset`, `denseStart/End`) and
`extractSeq` are in the **original read** (read `+` strand, 5'→3') frame.

### Type-1 — `C` (21 columns)

```
C  readName readLen extractStart extractEnd clipOffset extractSeq
   readStart1 readEnd1 strand1 ctg1 ctgStart1 ctgEnd1 mapq1
   readStart2 readEnd2 strand2 ctg2 ctgStart2 ctgEnd2 mapq2
```

- `extractSeq` is the read subsequence `[extractStart, extractEnd)` — up to `-l` bp
  either side of `clipOffset`.
- `clipOffset` is the middle of the clip:
  - **(a)** clipped sequence is not aligned elsewhere → the clip boundary;
  - **(b/c)** the neighbouring segment aligns (gap **or** overlap) →
    `(readEnd1 + readStart2) / 2`.
- alignment 1 = the segment *before* the clip on the read, alignment 2 = the
  segment *after*; segments are ordered by their position on the read. A missing
  neighbour (a terminal clip) is written as `-1` (integers) and `*` (strand/contig).

### Type-2 — `D` (15 columns)

```
D  readName readLen extractStart extractEnd denseStart denseEnd extractSeq
   readStart readEnd strand ctg ctgStart ctgEnd mapq
```

- `denseStart/denseEnd` bound the mismatch-dense region (first→last high-quality
  mismatch of a qualifying stretch; overlapping windows are merged).
- `extractStart/End` add `-l` bp of flanking either side of the dense region.
- The remaining columns describe the alignment the region belongs to.

## Type-2 requirements

A `D` line needs both base qualities **and** a source of mismatch positions. The
source is chosen in priority order:

1. reference genome via `-r`,
2. a `cs` tag,
3. an `MD` tag,
4. `X` (sequence-mismatch) ops in the CIGAR.

If base qualities are absent (e.g. many ONT BAMs) or no source is available, no
`D` lines are emitted (a single warning is printed to stderr for the latter).

## Tests

```sh
cargo test
```
