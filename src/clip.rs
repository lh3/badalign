//! Type-1 (C) lines: long clips at the primary alignment's clip boundaries.

use crate::align::{extract, Aln};
use std::io::{self, Write};

/// Emit up to two C-lines for a primary alignment.
///
/// `fwd_read` is the read in original (read `+`) orientation; `off` is the
/// read-forward coordinate of `fwd_read[0]` (0 unless the primary is hard-clipped).
#[allow(clippy::too_many_arguments)]
pub fn emit_c_lines<W: Write>(
    w: &mut W,
    name: &str,
    seg: u8,
    primary: &Aln,
    sa: &[Aln],
    fwd_read: &[u8],
    off: usize,
    read_len: usize,
    min_clip: usize,
    flank: usize,
) -> io::Result<()> {
    // All segments ordered along the read.
    let mut alns: Vec<&Aln> = Vec::with_capacity(1 + sa.len());
    alns.push(primary);
    alns.extend(sa.iter());
    alns.sort_by_key(|a| (a.q_start, a.q_end));
    let pi = alns.iter().position(|a| a.is_primary).unwrap();

    let clip5 = primary.q_start;
    let clip3 = read_len - primary.q_end;

    if clip5 >= min_clip {
        let prev = if pi > 0 { Some(alns[pi - 1]) } else { None };
        let clip_offset = match prev {
            None => primary.q_start,
            Some(p) => (p.q_end + primary.q_start) / 2,
        };
        emit_one(
            w,
            name,
            seg,
            read_len,
            fwd_read,
            off,
            flank,
            clip_offset,
            prev,
            Some(primary),
        )?;
    }
    if clip3 >= min_clip {
        let next = if pi + 1 < alns.len() {
            Some(alns[pi + 1])
        } else {
            None
        };
        let clip_offset = match next {
            None => primary.q_end,
            Some(n) => (primary.q_end + n.q_start) / 2,
        };
        emit_one(
            w,
            name,
            seg,
            read_len,
            fwd_read,
            off,
            flank,
            clip_offset,
            Some(primary),
            next,
        )?;
    }
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn emit_one<W: Write>(
    w: &mut W,
    name: &str,
    seg: u8,
    read_len: usize,
    fwd_read: &[u8],
    off: usize,
    flank: usize,
    clip_offset: usize,
    align1: Option<&Aln>,
    align2: Option<&Aln>,
) -> io::Result<()> {
    let want_start = clip_offset.saturating_sub(flank);
    let want_end = (clip_offset + flank).min(read_len);
    let (estart, eend, seq) = match extract(fwd_read, off, want_start, want_end) {
        Some(v) => v,
        None => return Ok(()),
    };
    // C<TAB>id<TAB>readName<TAB>... where id = readName_seg_C_extractStart_extractEnd.
    write!(
        w,
        "C\t{name}_{seg}_C_{estart}_{eend}\t{name}\t{read_len}\t{estart}\t{eend}\t{clip_offset}\t"
    )?;
    w.write_all(&seq)?;
    write!(w, "\t")?;
    write_aln(w, align1)?;
    write!(w, "\t")?;
    write_aln(w, align2)?;
    writeln!(w)?;
    Ok(())
}

/// Seven columns describing one alignment, or placeholders when absent.
fn write_aln<W: Write>(w: &mut W, a: Option<&Aln>) -> io::Result<()> {
    match a {
        Some(a) => write!(
            w,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}",
            a.q_start, a.q_end, a.strand as char, a.ctg, a.ref_start, a.ref_end, a.mapq
        ),
        None => write!(w, "-1\t-1\t*\t*\t-1\t-1\t-1"),
    }
}
