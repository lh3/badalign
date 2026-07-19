# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

`badalign` is a Rust CLI that scans a BAM and emits two kinds of "bad alignment" records (TAB-delimited, to stdout):

- **`C` lines (type-1)** — long clips (≥`-c`, default 40 bp) at the **primary** alignment's clip boundaries.
- **`D` lines (type-2)** — mismatch-dense regions (≥`-m`, default 10, high-quality mismatches within a `-w`=100 bp window).

Built on the `noodles` crate for BAM/FASTA I/O. See `README.md` for the exact output column layouts.

## Commands

```sh
cargo build --release          # binary at target/release/badalign
cargo test                     # unit tests (in-module #[cfg(test)] blocks)
cargo test dense::tests::single_region_first_to_last   # a single test by path
cargo fmt --check              # CI-style format check (run `cargo fmt` to fix)
cargo clippy --release         # lints; keep this at zero warnings
```

Run against the sample BAMs under `/Users/hli/work/minibwa/` (developer machine):
`HG002.HiFi-1k.bam` is smallest/fastest and has a `cs` tag (good for both C and D lines);
`HG002.HiFi-10k.bam` also has `cs`; `HG002.ONT-10k.bam` has **no base quality** (C lines only);
`hs38.fa` is the reference for the `-r` path. D-line changes should be validated against an
independent recomputation (see below), not just unit tests.

## Architecture — the load-bearing ideas

Processing is **streaming and order-agnostic** (`main.rs::process_record` is called per record; there is no per-QNAME buffering). It must keep working on **coordinate-sorted** BAMs, so never introduce logic that assumes a read's records are adjacent.

**The read-forward coordinate frame is the crux.** Every read-based number in the output (`readStart/End`, `extractStart/End`, `clipOffset`, `denseStart/End`, and `extractSeq`) is in the *original read* (read `+` strand, 5'→3') frame. This frame is derived **per-record** from `(strand, CIGAR clips, read length)` in `align.rs` (`compute_geom` → `query_interval`), so it never needs the primary record present. Two formulas encode it and must stay consistent:
- segment interval: `+` → `(lead_clip, L-tail_clip)`, `-` → `(tail_clip, L-lead_clip)`;
- SEQ-index → read-forward position for a mismatch: `+` → `lead_hard + rec_q`, `-` → `L-1-lead_hard-rec_q`.

`L` (read length) counts hard clips (`M+I+S+=+X+H`); `rec_q` indexes the record SEQ (which also indexes QUAL). `extract()` reconstructs subsequences by clamping to the bases actually available in a (possibly hard-clipped) record, using an `off` = leading (`+`) / trailing (`-`) hard-clip length.

**C lines are generated only while processing the primary record** (`clip.rs`); secondary/supplementary records are skipped. The segment set is the primary plus the segments reconstructed from the primary's **`SA` tag** (`align.rs::parse_sa`), *not* the supplementary records — deliberate so it works regardless of sort order. Segments shorter than `-f`×`readLen` on the read are ignored first, then `clip::remove_contained` drops any segment whose read interval is strictly contained in a longer one (both uniform — the primary itself may be dropped). A C line is then emitted at each boundary of the remaining read-sorted chain: the two **terminal** clips (gated by `-c` min_clip) and each **internal** junction between consecutive segments (gated by `-g` max gap/overlap). `emit_one` builds the outer-span extract, the `leftflank/rightflank` id fields, and writes the line; a terminal clip's absent neighbour prints `-1`/`*` placeholders.

On a D line, `ctgStart/ctgEnd` are the **dense region's** reference span (min→max event ref, half-open), not the whole alignment — mismatch events carry `ref_pos` (`mismatch.rs`) which `hq_positions` keeps as `(fwd_pos, ref_pos)` pairs so `emit_d_lines` can derive it.

**D lines come from every mapped record** with base quality (`dense.rs`). Substitution positions come from a source chosen in priority order in `main.rs::choose_sites`: reference (`-r`) > `cs` > `MD` > CIGAR `X` (`mismatch.rs::from_reference`/`from_cs`/`from_md`/`from_cigar_x`). **Gate:** a D line requires base quality AND one of those substitution sources — do not loosen this (e.g. gap-only D lines were explicitly declined). Gap opens (`mismatch.rs::gap_events`) are added to the count only for records that already pass the gate; they come from the CIGAR (`I`/`D`, one event each, `N` excluded) with quality = max over the gap plus flanking bases. The sliding-window merge that turns events into regions is `dense.rs::dense_regions`.

## Conventions / invariants to preserve

- All emitted positions are **0-based, half-open** (BED convention).
- Every line starts with a type letter (`C`/`D`) then an **id** built from its own fields (`<seg>` = 0 single-end / 1,2 mates): C = `readName_<seg>_C_<extractStart>_<extractEnd>_<leftflank>_<rightflank>` (leftflank = readEnd1−extractStart, rightflank = extractEnd−readStart2; 0 if that neighbour is absent, signed), D = `readName_<seg>_D_<extractStart>_<extractEnd>`. Because the id is field-derived, the emit functions write directly and output streams in BAM order — no buffering. It may rarely collide (two lines of a read+mate+type with the same extract window).
- Column counts are fixed: **C = 22**, **D = 17**. Adding/removing a field breaks downstream parsing — update `README.md` and the output writers together.
- When touching D-line counting, verify by independently recomputing events (substitutions + gap opens with the max-flank quality rule) and the window/merge, comparing `denseStart/denseEnd` per record across both strands and hard-clipped supplementaries — unit tests alone won't catch coordinate-frame regressions.
- CIGAR is handled as a flat `Vec<(usize, u8)>` (`align::CigarOps`) with op letters `M I D N S H P = X`; `main::kind_byte` maps noodles' `Kind` into these.
