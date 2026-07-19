//! Type-1 (C) lines: junctions of a read's alignment chain.
//!
//! Emitted only while processing the **primary** record: the segment set is the
//! primary plus the segments parsed from its `SA` tag (supplementary records are
//! never read directly). Contained segments are dropped, then a C-line is emitted
//! at each boundary of the remaining, read-sorted chain — the two terminal clips
//! (gated by `-c`) and each internal junction (gated by `-g`).

use crate::align::{extract, Aln};
use std::io::{self, Write};

/// Drop any segment strictly contained (on read/query coords) in a longer one.
fn remove_contained<'a>(alns: &[&'a Aln]) -> Vec<&'a Aln> {
    let n = alns.len();
    (0..n)
        .filter(|&i| {
            !(0..n).any(|j| {
                j != i
                    && alns[j].q_start <= alns[i].q_start
                    && alns[i].q_end <= alns[j].q_end
                    && (alns[j].q_start < alns[i].q_start || alns[j].q_end > alns[i].q_end)
            })
        })
        .map(|i| alns[i])
        .collect()
}

/// Emit C-lines for a read, given its primary alignment and `SA` segments.
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
    max_gap: usize,
    min_frac: f64,
    flank: usize,
) -> io::Result<()> {
    let mut all: Vec<&Aln> = Vec::with_capacity(1 + sa.len());
    all.push(primary);
    all.extend(sa.iter());
    // Ignore alignments whose read span is shorter than min_frac * readLen.
    let min_len = min_frac * read_len as f64;
    all.retain(|a| (a.q_end - a.q_start) as f64 >= min_len);
    // Filter out contained alignments (uniformly — the primary may be removed),
    // then order the remaining chain along the read.
    let mut chain = remove_contained(&all);
    if chain.is_empty() {
        return Ok(());
    }
    chain.sort_by_key(|a| (a.q_start, a.q_end));

    // 5' terminal clip.
    let first = chain[0];
    if first.q_start >= min_clip {
        emit_one(
            w,
            name,
            seg,
            read_len,
            fwd_read,
            off,
            flank,
            first.q_start,
            None,
            Some(first),
        )?;
    }
    // Internal junctions (gated by max gap/overlap).
    for pair in chain.windows(2) {
        let (a, b) = (pair[0], pair[1]);
        let gap_overlap = (b.q_start as i64 - a.q_end as i64).unsigned_abs() as usize;
        if gap_overlap <= max_gap {
            let clip_offset = (a.q_end + b.q_start) / 2;
            emit_one(
                w,
                name,
                seg,
                read_len,
                fwd_read,
                off,
                flank,
                clip_offset,
                Some(a),
                Some(b),
            )?;
        }
    }
    // 3' terminal clip.
    let last = chain[chain.len() - 1];
    if read_len - last.q_end >= min_clip {
        emit_one(
            w,
            name,
            seg,
            read_len,
            fwd_read,
            off,
            flank,
            last.q_end,
            Some(last),
            None,
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
    // The extract spans the whole junction: `flank` bp outside the previous
    // alignment's end and the next alignment's start. Using min/max of the two
    // bounds covers the entire gap (readStart2 > readEnd1) or overlap
    // (readStart2 < readEnd1) region without ever inverting. For a terminal clip
    // (a neighbour is absent) the bound falls back to the clip boundary
    // (clip_offset), giving a `2*flank` window around it.
    let left_bound = align1.map_or(clip_offset, |a| a.q_end);
    let right_bound = align2.map_or(clip_offset, |a| a.q_start);
    let want_start = left_bound.min(right_bound).saturating_sub(flank);
    let want_end = (left_bound.max(right_bound) + flank).min(read_len);
    let (estart, eend, seq) = match extract(fwd_read, off, want_start, want_end) {
        Some(v) => v,
        None => return Ok(()),
    };
    // id = readName_seg_C_extractStart_extractEnd_leftflank_rightflank, where
    // leftflank/rightflank measure how much of the extract each neighbour covers
    // (0 when that neighbour is absent; signed — a large gap can make them negative).
    let leftflank = align1.map_or(0, |a| a.q_end as i64 - estart as i64);
    let rightflank = align2.map_or(0, |a| eend as i64 - a.q_start as i64);
    write!(
        w,
        "C\t{name}_{seg}_C_{estart}_{eend}_{leftflank}_{rightflank}\t{name}\t{read_len}\t{estart}\t{eend}\t{clip_offset}\t"
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

#[cfg(test)]
mod tests {
    use super::*;

    fn aln(q_start: usize, q_end: usize) -> Aln {
        Aln {
            ctg: "c".to_string(),
            ref_start: 0,
            ref_end: q_end - q_start,
            strand: b'+',
            mapq: 60,
            q_start,
            q_end,
        }
    }

    #[test]
    fn contained_segment_removed() {
        let a = aln(0, 1000); // container
        let b = aln(200, 300); // contained in a
        let c = aln(1000, 2000); // adjacent, not contained
        let all: Vec<&Aln> = vec![&a, &b, &c];
        let kept: Vec<(usize, usize)> = remove_contained(&all)
            .iter()
            .map(|x| (x.q_start, x.q_end))
            .collect();
        assert_eq!(kept, vec![(0, 1000), (1000, 2000)]); // b dropped, order preserved
    }

    #[test]
    fn identical_intervals_both_kept() {
        let a = aln(0, 500);
        let b = aln(0, 500); // same span (e.g. a multi-mapping) -> neither strictly contains
        let all: Vec<&Aln> = vec![&a, &b];
        assert_eq!(remove_contained(&all).len(), 2);
    }

    #[test]
    fn nested_containment() {
        let a = aln(0, 1000);
        let b = aln(100, 900); // in a
        let c = aln(200, 800); // in a and b
        let all: Vec<&Aln> = vec![&a, &b, &c];
        let kept: Vec<(usize, usize)> = remove_contained(&all)
            .iter()
            .map(|x| (x.q_start, x.q_end))
            .collect();
        assert_eq!(kept, vec![(0, 1000)]);
    }
}
