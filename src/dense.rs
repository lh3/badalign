//! Type-2 (D) lines: mismatch-dense regions of one alignment record.

use crate::align::{extract, Aln};
use std::io::{self, Write};

/// Map mismatch events `(rec_q, baseq)` to sorted read-forward positions, keeping
/// only high-quality ones. `rec_q` indexes the record SEQ; `lead_hard` is the
/// record's leading hard-clip length (0 for a primary).
pub fn hq_positions(
    events: &[(usize, u8)],
    strand: u8,
    read_len: usize,
    lead_hard: usize,
    min_baseq: u8,
) -> Vec<usize> {
    let mut pos: Vec<usize> = events
        .iter()
        .filter_map(|&(rec_q, bq)| {
            if bq < min_baseq {
                return None;
            }
            let fwd = if strand == b'+' {
                lead_hard + rec_q
            } else {
                read_len - 1 - lead_hard - rec_q
            };
            Some(fwd)
        })
        .collect();
    pos.sort_unstable();
    pos
}

/// Maximal regions where some `w`-wide window holds `>= m` mismatches. Each
/// region spans from its first to its last mismatch (`denseEnd` is exclusive).
pub fn dense_regions(positions: &[usize], w: usize, m: usize) -> Vec<(usize, usize)> {
    let mut regions = Vec::new();
    if positions.is_empty() {
        return regions;
    }
    let mut front = 0usize;
    let mut reg: Option<(usize, usize)> = None; // (start, end_inclusive)
    for i in 0..positions.len() {
        let p = positions[i];
        while positions[front] + w <= p {
            front += 1;
        }
        if i - front + 1 >= m {
            let wstart = positions[front];
            reg = Some(match reg {
                Some((s, e)) if wstart <= e => (s, e.max(p)),
                Some((s, e)) => {
                    regions.push((s, e + 1));
                    (wstart, p)
                }
                None => (wstart, p),
            });
        }
    }
    if let Some((s, e)) = reg {
        regions.push((s, e + 1));
    }
    regions
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn no_region_below_threshold() {
        // 4 mismatches in a 100bp window, m=5 -> nothing.
        let pos = vec![10, 30, 50, 70];
        assert!(dense_regions(&pos, 100, 5).is_empty());
    }

    #[test]
    fn single_region_first_to_last() {
        let pos = vec![10, 20, 30, 40, 50];
        assert_eq!(dense_regions(&pos, 100, 5), vec![(10, 51)]);
    }

    #[test]
    fn continuing_windows_merge() {
        let pos = vec![10, 20, 30, 40, 50, 60];
        assert_eq!(dense_regions(&pos, 100, 5), vec![(10, 61)]);
    }

    #[test]
    fn separate_clusters_split() {
        // Two dense clusters far apart -> two regions.
        let mut pos = vec![0, 10, 20, 30, 40];
        pos.extend([500, 510, 520, 530, 540]);
        assert_eq!(dense_regions(&pos, 100, 5), vec![(0, 41), (500, 541)]);
    }

    #[test]
    fn window_width_excludes_far_mismatch() {
        // Spread over >100bp: no 100bp window holds 5.
        let pos = vec![0, 30, 60, 90, 120];
        assert!(dense_regions(&pos, 100, 5).is_empty());
    }
}

/// Emit D-lines for a record given its high-quality mismatch positions.
#[allow(clippy::too_many_arguments)]
pub fn emit_d_lines<W: Write>(
    w: &mut W,
    name: &str,
    aln: &Aln,
    positions: &[usize],
    fwd_seq: &[u8],
    off: usize,
    read_len: usize,
    window: usize,
    min_mismatch: usize,
    flank: usize,
) -> io::Result<()> {
    for (dstart, dend) in dense_regions(positions, window, min_mismatch) {
        let want_start = dstart.saturating_sub(flank);
        let want_end = (dend + flank).min(read_len);
        let (estart, eend, seq) = match extract(fwd_seq, off, want_start, want_end) {
            Some(v) => v,
            None => continue,
        };
        // Count of mismatches + gap opens inside the dense interval (positions is
        // sorted; flanking events are excluded).
        let n =
            positions.partition_point(|&p| p < dend) - positions.partition_point(|&p| p < dstart);
        write!(
            w,
            "D\t{name}\t{read_len}\t{estart}\t{eend}\t{dstart}\t{dend}\t{n}\t"
        )?;
        w.write_all(&seq)?;
        writeln!(
            w,
            "\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            aln.q_start,
            aln.q_end,
            aln.strand as char,
            aln.ctg,
            aln.ref_start,
            aln.ref_end,
            aln.mapq
        )?;
    }
    Ok(())
}
