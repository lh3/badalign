mod align;
mod cli;
mod clip;
mod dense;
mod mismatch;
mod refgenome;

use align::{compute_geom, parse_sa, query_interval, revcomp, Aln};
use anyhow::{Context, Result};
use clap::Parser;
use noodles::sam::alignment::record::cigar::op::Kind;
use noodles::sam::alignment::record::data::field::Value;
use refgenome::RefGenome;
use std::io::{self, BufWriter, Write};
use std::sync::Once;

fn main() -> Result<()> {
    let args = cli::Args::parse();

    // The reference is only a D-line mismatch source, so skip loading it (and its
    // ~3 GB of RAM) unless D-lines are requested.
    let refgenome = match (&args.reference, args.emit_d) {
        (Some(p), true) => Some(RefGenome::load(p).context("loading reference genome")?),
        _ => None,
    };

    let mut reader = noodles::bam::io::reader::Builder
        .build_from_path(&args.bam)
        .with_context(|| format!("opening BAM {}", args.bam.display()))?;
    let header = reader.read_header().context("reading BAM header")?;
    let ref_names: Vec<String> = header
        .reference_sequences()
        .keys()
        .map(|k| String::from_utf8_lossy(k).into_owned())
        .collect();

    let stdout = io::stdout();
    let mut out = BufWriter::new(stdout.lock());
    for result in reader.records() {
        let record = result.context("reading BAM record")?;
        process_record(&mut out, &record, &ref_names, refgenome.as_ref(), &args)?;
    }
    out.flush()?;
    Ok(())
}

fn kind_byte(k: Kind) -> u8 {
    match k {
        Kind::Match => b'M',
        Kind::Insertion => b'I',
        Kind::Deletion => b'D',
        Kind::Skip => b'N',
        Kind::SoftClip => b'S',
        Kind::HardClip => b'H',
        Kind::Pad => b'P',
        Kind::SequenceMatch => b'=',
        Kind::SequenceMismatch => b'X',
    }
}

fn get_str_tag(record: &noodles::bam::Record, tag: &[u8; 2]) -> Option<Vec<u8>> {
    match record.data().get(tag) {
        Some(Ok(Value::String(s))) => Some(s.to_vec()),
        _ => None,
    }
}

fn warn_no_source() {
    static ONCE: Once = Once::new();
    ONCE.call_once(|| {
        eprintln!(
            "badalign: no mismatch source (need -r, or a cs/MD tag, or X in CIGAR); \
             skipping D-lines for such records"
        );
    });
}

fn process_record<W: Write>(
    out: &mut W,
    record: &noodles::bam::Record,
    ref_names: &[String],
    refgenome: Option<&RefGenome>,
    args: &cli::Args,
) -> Result<()> {
    let flags = record.flags();
    if flags.is_unmapped() || flags.is_secondary() {
        return Ok(());
    }
    let is_primary = !flags.is_supplementary();

    // Segment: 0 = single-end, 1/2 = first/second mate of a pair.
    let seg: u8 = if !flags.is_segmented() {
        0
    } else if flags.is_first_segment() {
        1
    } else if flags.is_last_segment() {
        2
    } else {
        0
    };

    // CIGAR -> flat ops.
    let mut ops = Vec::new();
    for op in record.cigar().iter() {
        let op = op.context("decoding CIGAR")?;
        ops.push((op.len(), kind_byte(op.kind())));
    }
    if ops.is_empty() {
        return Ok(());
    }
    let g = compute_geom(&ops);
    let read_len = g.read_len;

    let strand = if flags.is_reverse_complemented() {
        b'-'
    } else {
        b'+'
    };
    let ref_id = record
        .reference_sequence_id()
        .context("missing reference id")??;
    let ctg = ref_names
        .get(ref_id)
        .cloned()
        .unwrap_or_else(|| ref_id.to_string());
    let ref_start = record
        .alignment_start()
        .context("missing alignment start")??
        .get()
        - 1;
    let ref_end = ref_start + g.ref_span;
    let mapq = record
        .mapping_quality()
        .map(|m| m.get() as i32)
        .unwrap_or(-1);
    let (q_start, q_end) = query_interval(&g, strand, read_len);

    let aln = Aln {
        ctg,
        ref_start,
        ref_end,
        strand,
        mapq,
        q_start,
        q_end,
    };

    let name = record
        .name()
        .map(|b| String::from_utf8_lossy(b).into_owned())
        .unwrap_or_else(|| "*".to_string());

    // Decoded record SEQ (ref-forward) and read-forward orientation.
    let seq: Vec<u8> = record.sequence().iter().collect();
    let fwd_seq = if strand == b'+' {
        seq.clone()
    } else {
        revcomp(&seq)
    };
    let off = if strand == b'+' {
        g.lead_hard
    } else {
        g.tail_hard
    };

    // --- Type-1 (C) lines: primary only ---
    if is_primary && !fwd_seq.is_empty() {
        let sa = get_str_tag(record, b"SA")
            .map(|s| parse_sa(&s, read_len))
            .unwrap_or_default();
        clip::emit_c_lines(
            out,
            &name,
            seg,
            &aln,
            &sa,
            &fwd_seq,
            off,
            read_len,
            args.min_clip,
            args.max_gap,
            args.min_frac,
            args.flank,
        )?;
    }

    // --- Type-2 (D) lines: opt-in via -d; any record with base quality + a source ---
    let qual = record.quality_scores();
    let qual = qual.as_bytes();
    if args.emit_d && !qual.is_empty() {
        let sites = choose_sites(record, &ops, &seq, ref_start, &aln.ctg, refgenome, &g);
        match sites {
            Some(sites) => {
                // Combine substitutions (quality at their base) with gap opens
                // (quality = max over the gap + flanking bases). Each event carries
                // its reference position so the dense region's ref span can be found.
                let mut events: Vec<(usize, usize, u8)> = sites
                    .iter()
                    .map(|s| (s.rec_q, s.ref_pos, qual.get(s.rec_q).copied().unwrap_or(0)))
                    .collect();
                events.extend(mismatch::gap_events(&ops, ref_start, qual));
                let positions =
                    dense::hq_positions(&events, strand, read_len, g.lead_hard, args.min_baseq);
                dense::emit_d_lines(
                    out,
                    &name,
                    seg,
                    &aln,
                    &positions,
                    &fwd_seq,
                    off,
                    read_len,
                    args.window,
                    args.min_mismatch,
                    args.flank,
                )?;
            }
            None => warn_no_source(),
        }
    }

    Ok(())
}

/// Pick a mismatch source in priority order: reference > cs > MD > CIGAR X.
fn choose_sites(
    record: &noodles::bam::Record,
    ops: &align::CigarOps,
    seq: &[u8],
    ref_start: usize,
    ctg: &str,
    refgenome: Option<&RefGenome>,
    g: &align::Geom,
) -> Option<Vec<mismatch::Site>> {
    if let Some(rg) = refgenome {
        if let Some(refseq) = rg.get(ctg) {
            return Some(mismatch::from_reference(ops, seq, ref_start, refseq));
        }
    }
    if let Some(cs) = get_str_tag(record, b"cs") {
        return Some(mismatch::from_cs(&cs, g.lead_soft, ref_start));
    }
    if let Some(md) = get_str_tag(record, b"MD") {
        return Some(mismatch::from_md(&md, ops, ref_start));
    }
    if g.has_x {
        return Some(mismatch::from_cigar_x(ops));
    }
    None
}
