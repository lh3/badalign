> [!Warning]
> This project is vibe coded with Claude Code.

# badalign

Extract *bad* read alignments from a BAM file. Two kinds of problems are reported,
both as TAB-delimited lines:

- **Type-1 (`C`)** — junctions of a read's alignment chain (long clips and chimeric
  breakpoints). Useful for finding chimeric/misassembled reads and breakpoints.
- **Type-2 (`D`)** — mismatch-dense regions (default ≥10 high-quality ≥Q20
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
  -g INT    max gap or overlap between adjacent alignments for a C-line (type-1) [100]
  -f FLOAT  ignore an alignment shorter than FLOAT*readLen on the read (type-1) [0.2]
  -m INT    min mismatches in a window (type-2) [10]
```

Output goes to stdout.

The input need not be sorted in any particular way — records are processed
streaming, and a read's chimeric segments are reconstructed from the primary
record's `SA` tag, so coordinate-sorted BAMs work as well as query-grouped ones.

## Output format

All coordinates are **0-based, half-open** (BED convention). Read-based positions
(`readStart/End`, `extractStart/End`, `clipOffset`, `denseStart/End`) and
`extractSeq` are in the **original read** (read `+` strand, 5'→3') frame.

Every line's first data column (before `readName`) is an **id**:
- C: `readName_<seg>_C_<extractStart>_<extractEnd>_<leftflank>_<rightflank>`
- D: `readName_<seg>_D_<extractStart>_<extractEnd>`

where `<seg>` is `0` single-end / `1`/`2` first/second mate of a paired-end read. For C,
`leftflank = readEnd1 − extractStart` and `rightflank = extractEnd − readStart2` measure how
much of the extract the previous/next alignment covers (`0` when that neighbour is absent, and
signed — a large gap can make a flank negative). `leftflank + rightflank` exceeds the extract
length exactly when the two adjacent alignments overlap on the read.

The id is self-contained, so output streams in BAM order (no buffering). It may collide
only if two lines of the same read+mate+type share the same extract window (very rare).

### Type-1 — `C` (22 columns)

```
C  id readName readLen extractStart extractEnd clipOffset extractSeq
   readStart1 readEnd1 strand1 ctg1 ctgStart1 ctgEnd1 mapq1
   readStart2 readEnd2 strand2 ctg2 ctgStart2 ctgEnd2 mapq2
```

- `extractSeq` is the read subsequence `[extractStart, extractEnd)`, spanning the whole
  junction: `-l` bp outside `readEnd1` and `readStart2`, i.e.
  `[min(readEnd1,readStart2) − l, max(readEnd1,readStart2) + l)`. This covers the entire
  unaligned gap (or the overlapping region) plus `-l` of aligned flank on each side. For a
  terminal clip (no neighbour on one side) it is `±-l` bp around the clip boundary.
- `clipOffset` is the middle of the clip:
  - **(a)** clipped sequence is not aligned elsewhere → the clip boundary;
  - **(b/c)** the neighbouring segment aligns (gap **or** overlap) →
    `(readEnd1 + readStart2) / 2`.
- alignment 1 = the segment *before* the clip on the read, alignment 2 = the
  segment *after*; segments are ordered by their position on the read. A missing
  neighbour (a terminal clip) is written as `-1` (integers) and `*` (strand/contig).

The read's segments are the primary alignment plus the segments in its `SA` tag
(C-lines are produced only while processing the primary record; supplementary/secondary
records are not read directly). A segment shorter than `-f`×`readLen` on the read is
ignored, and a segment whose read interval is contained in a longer one is dropped. A C-line is then emitted at each boundary of the remaining segment chain:
the two **terminal** clips (each reported when its clip length is ≥ `-c`) and each
**internal** junction between consecutive segments (reported when the gap or overlap
between them is ≤ `-g` — this suppresses junctions whose extract would span most of the
read).

### Type-2 — `D` (17 columns)

```
D  id readName readLen extractStart extractEnd denseStart denseEnd nMismatch extractSeq
   readStart readEnd strand ctg ctgStart ctgEnd mapq
```

- `denseStart/denseEnd` bound the mismatch-dense region **in read coordinates**
  (first→last high-quality mismatch of a qualifying stretch; overlapping windows merged).
- `ctgStart/ctgEnd` are the same dense region **in reference coordinates** (min→max
  reference position of its events, half-open) — not the whole alignment's span.
- `nMismatch` is the number of mismatches + gap opens inside the dense interval
  `[denseStart, denseEnd)` (flanking events are not counted).
- `extractStart/End` add `-l` bp of flanking either side of the dense region.
- `readStart/readEnd`, `strand`, `ctg`, `mapq` describe the alignment the region belongs to.

A "mismatch" here is either a **substitution** (base quality ≥ `-Q` at that base) or
a **gap open** — one event per insertion/deletion regardless of length (gap
*extensions* don't count; introns/`N` are not gaps). A gap's base quality is the
highest quality among the bases in the gap plus the read bases immediately before
and after it (a deletion has no bases in the gap, so just its two flanks), and it
counts when that quality is ≥ `-Q`.

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
