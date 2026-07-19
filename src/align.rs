//! Alignment geometry: CIGAR parsing, read-forward coordinate math, and SA-tag
//! reconstruction of a read's chimeric segments.
//!
//! All read coordinates are expressed in the *original read* (read `+` strand,
//! 5'->3') frame, derived per-record from (strand, CIGAR clips, read length) so
//! it never needs the primary record to be present.

/// A CIGAR as a flat list of (length, op-char) where op-char is one of
/// `M I D N S H P = X` (the SAM CIGAR letters).
pub type CigarOps = Vec<(usize, u8)>;

/// One alignment segment of a read (primary or an `SA` segment).
#[derive(Clone, Debug)]
pub struct Aln {
    pub ctg: String,
    pub ref_start: usize, // 0-based, half-open
    pub ref_end: usize,
    pub strand: u8,     // b'+' or b'-'
    pub mapq: i32,      // -1 if missing
    pub q_start: usize, // read-forward coords, 0-based half-open
    pub q_end: usize,
}

/// Geometry derived from a single CIGAR.
pub struct Geom {
    pub lead_clip: usize, // leading S+H (read bases before the aligned block)
    pub tail_clip: usize, // trailing S+H
    pub lead_hard: usize, // leading H only
    pub tail_hard: usize, // trailing H only
    pub lead_soft: usize, // leading S only
    pub ref_span: usize,  // reference bases consumed
    #[allow(dead_code)]
    pub seq_len: usize, // bases stored in SEQ (M+I+S+=+X)
    pub read_len: usize,  // full read length (seq_len + lead_hard + tail_hard)
    pub has_x: bool,      // contains a SequenceMismatch (X) op
}

pub fn compute_geom(ops: &CigarOps) -> Geom {
    let mut lead_clip = 0;
    let mut lead_hard = 0;
    let mut lead_soft = 0;
    // leading clip run
    for &(len, k) in ops {
        match k {
            b'H' => {
                lead_clip += len;
                lead_hard += len;
            }
            b'S' => {
                lead_clip += len;
                lead_soft += len;
            }
            _ => break,
        }
    }
    let mut tail_clip = 0;
    let mut tail_hard = 0;
    for &(len, k) in ops.iter().rev() {
        match k {
            b'H' => {
                tail_clip += len;
                tail_hard += len;
            }
            b'S' => {
                tail_clip += len;
            }
            _ => break,
        }
    }
    let mut ref_span = 0;
    let mut seq_len = 0;
    let mut has_x = false;
    for &(len, k) in ops {
        match k {
            b'M' | b'=' | b'X' => {
                ref_span += len;
                seq_len += len;
                if k == b'X' {
                    has_x = true;
                }
            }
            b'I' | b'S' => seq_len += len,
            b'D' | b'N' => ref_span += len,
            _ => {} // H, P
        }
    }
    let read_len = seq_len + lead_hard + tail_hard;
    Geom {
        lead_clip,
        tail_clip,
        lead_hard,
        tail_hard,
        lead_soft,
        ref_span,
        seq_len,
        read_len,
        has_x,
    }
}

/// Read-forward query interval `[q_start, q_end)` for a segment, given its
/// leading/trailing clips, strand, and the full read length.
pub fn query_interval(g: &Geom, strand: u8, read_len: usize) -> (usize, usize) {
    if strand == b'+' {
        (g.lead_clip, read_len - g.tail_clip)
    } else {
        (g.tail_clip, read_len - g.lead_clip)
    }
}

/// Parse a SAM CIGAR string (as found in an `SA` tag) into ops.
pub fn parse_cigar_str(s: &[u8]) -> Option<CigarOps> {
    let mut ops = Vec::new();
    let mut n: usize = 0;
    let mut seen = false;
    for &b in s {
        if b.is_ascii_digit() {
            n = n * 10 + (b - b'0') as usize;
            seen = true;
        } else {
            if !seen {
                return None;
            }
            ops.push((n, b));
            n = 0;
            seen = false;
        }
    }
    if seen {
        return None;
    }
    Some(ops)
}

/// Parse an `SA:Z` value into segments, given the read length taken from the
/// primary. Format: `rname,pos,strand,CIGAR,mapQ,NM;` repeated.
pub fn parse_sa(sa: &[u8], read_len: usize) -> Vec<Aln> {
    let mut out = Vec::new();
    for seg in sa.split(|&b| b == b';') {
        if seg.is_empty() {
            continue;
        }
        let f: Vec<&[u8]> = seg.split(|&b| b == b',').collect();
        if f.len() < 5 {
            continue;
        }
        let ctg = String::from_utf8_lossy(f[0]).into_owned();
        let pos: usize = match std::str::from_utf8(f[1]).ok().and_then(|s| s.parse().ok()) {
            Some(p) => p,
            None => continue,
        };
        let strand = if f[2] == b"-" { b'-' } else { b'+' };
        let ops = match parse_cigar_str(f[3]) {
            Some(o) => o,
            None => continue,
        };
        let mapq: i32 = std::str::from_utf8(f[4])
            .ok()
            .and_then(|s| s.parse().ok())
            .unwrap_or(-1);
        let g = compute_geom(&ops);
        let ref_start = pos.saturating_sub(1);
        let ref_end = ref_start + g.ref_span;
        let (q_start, q_end) = query_interval(&g, strand, read_len);
        out.push(Aln {
            ctg,
            ref_start,
            ref_end,
            strand,
            mapq,
            q_start,
            q_end,
        });
    }
    out
}

/// Extract read-forward bases for `[start, end)` from `fwd_seq`, whose first base
/// sits at read-forward coordinate `off`. Returns the clamped `(start, end, seq)`
/// intersected with the bases actually available in `fwd_seq`.
pub fn extract(
    fwd_seq: &[u8],
    off: usize,
    start: usize,
    end: usize,
) -> Option<(usize, usize, Vec<u8>)> {
    let s = start.max(off);
    let e = end.min(off + fwd_seq.len());
    if s >= e {
        return None;
    }
    Some((s, e, fwd_seq[s - off..e - off].to_vec()))
}

pub fn revcomp(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&b| complement(b)).collect()
}

pub fn complement(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        b'a' => b't',
        b'c' => b'g',
        b'g' => b'c',
        b't' => b'a',
        _ => b'N',
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn geom_and_query_interval_forward() {
        // 10254S 100M 2S on a 10356 bp read, + strand.
        let ops = parse_cigar_str(b"10254S100M2S").unwrap();
        let g = compute_geom(&ops);
        assert_eq!(g.lead_clip, 10254);
        assert_eq!(g.tail_clip, 2);
        assert_eq!(g.read_len, 10254 + 100 + 2);
        assert_eq!(g.ref_span, 100);
        let (qs, qe) = query_interval(&g, b'+', g.read_len);
        assert_eq!((qs, qe), (10254, 10354));
    }

    #[test]
    fn query_interval_reverse_is_flipped() {
        // Same clips, but on '-' strand the leading/trailing swap in read space.
        let ops = parse_cigar_str(b"10S100M20S").unwrap();
        let g = compute_geom(&ops);
        let l = g.read_len; // 130
        let (qs, qe) = query_interval(&g, b'-', l);
        // '-': q_start = tail_clip = 20, q_end = L - lead_clip = 130 - 10 = 120
        assert_eq!((qs, qe), (20, 120));
    }

    #[test]
    fn hardclip_counts_toward_read_len() {
        // Supplementary style: 33H 49M 8H  => read_len 90, lead_hard 33.
        let ops = parse_cigar_str(b"33H49M8H").unwrap();
        let g = compute_geom(&ops);
        assert_eq!(g.lead_hard, 33);
        assert_eq!(g.tail_hard, 8);
        assert_eq!(g.read_len, 90);
        assert_eq!(g.seq_len, 49);
    }

    #[test]
    fn parse_sa_segment() {
        let read_len = 16920;
        let alns = parse_sa(b"chr17,25068925,+,33S8120M497D8767S,4,3464;", read_len);
        assert_eq!(alns.len(), 1);
        let a = &alns[0];
        assert_eq!(a.ctg, "chr17");
        assert_eq!(a.ref_start, 25068924);
        assert_eq!(a.ref_end, 25068924 + 8120 + 497);
        assert_eq!(a.strand, b'+');
        assert_eq!(a.mapq, 4);
        assert_eq!((a.q_start, a.q_end), (33, 16920 - 8767));
    }

    #[test]
    fn extract_clamps_to_available() {
        let fwd = b"ACGTACGTAC".to_vec();
        // off=5 means fwd[0] is read position 5; available window [5,15).
        let (s, e, seq) = extract(&fwd, 5, 3, 12).unwrap();
        assert_eq!((s, e), (5, 12));
        assert_eq!(seq, b"ACGTACG");
        assert!(extract(&fwd, 5, 0, 5).is_none());
    }

    #[test]
    fn revcomp_works() {
        assert_eq!(revcomp(b"ACGTN"), b"NACGT");
    }
}
