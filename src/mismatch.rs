//! Enumerate substitution sites of one alignment record from whichever source
//! is available, in priority order: reference FASTA > cs tag > MD tag > CIGAR X.
//!
//! A site is `(rec_q, ref_pos)` where `rec_q` indexes the record SEQ (ref-forward
//! orientation, so it also indexes the record's QUAL) and `ref_pos` is 0-based on
//! the contig. Indels are ignored — only substitutions count.

use crate::align::CigarOps;

pub struct Site {
    pub rec_q: usize,
    /// Reference position of the substitution (retained for future reference-based
    /// reporting; not currently emitted).
    #[allow(dead_code)]
    pub ref_pos: usize,
}

fn base_eq(a: u8, b: u8) -> bool {
    a.eq_ignore_ascii_case(&b)
}

/// Compare read vs reference over every query+ref-consuming op.
pub fn from_reference(ops: &CigarOps, seq: &[u8], ref_start: usize, refseq: &[u8]) -> Vec<Site> {
    let mut out = Vec::new();
    let mut q = 0usize;
    let mut r = ref_start;
    for &(len, k) in ops {
        match k {
            b'M' | b'=' | b'X' => {
                for i in 0..len {
                    let sb = seq.get(q + i).copied().unwrap_or(b'N');
                    if let Some(&rb) = refseq.get(r + i) {
                        if !rb.eq_ignore_ascii_case(&b'N')
                            && !sb.eq_ignore_ascii_case(&b'N')
                            && !base_eq(sb, rb)
                        {
                            out.push(Site {
                                rec_q: q + i,
                                ref_pos: r + i,
                            });
                        }
                    }
                }
                q += len;
                r += len;
            }
            b'I' | b'S' => q += len,
            b'D' | b'N' => r += len,
            _ => {} // H, P
        }
    }
    out
}

/// Mismatches straight from CIGAR `X` ops.
pub fn from_cigar_x(ops: &CigarOps) -> Vec<Site> {
    let mut out = Vec::new();
    let mut q = 0usize;
    let mut r = 0usize;
    for &(len, k) in ops {
        match k {
            b'X' => {
                for i in 0..len {
                    out.push(Site {
                        rec_q: q + i,
                        ref_pos: r + i,
                    });
                }
                q += len;
                r += len;
            }
            b'M' | b'=' => {
                q += len;
                r += len;
            }
            b'I' | b'S' => q += len,
            b'D' | b'N' => r += len,
            _ => {}
        }
    }
    out
}

/// Parse a `cs:Z` string (short or long form). `lead_soft` is the leading soft
/// clip (SEQ offset of the first aligned base); `ref_start` is the alignment
/// start on the contig.
pub fn from_cs(cs: &[u8], lead_soft: usize, ref_start: usize) -> Vec<Site> {
    let mut out = Vec::new();
    let mut q = lead_soft;
    let mut r = ref_start;
    let mut i = 0;
    while i < cs.len() {
        let op = cs[i];
        i += 1;
        match op {
            b':' => {
                // :N  -> N identical bases
                let mut n = 0usize;
                while i < cs.len() && cs[i].is_ascii_digit() {
                    n = n * 10 + (cs[i] - b'0') as usize;
                    i += 1;
                }
                q += n;
                r += n;
            }
            b'*' => {
                // *xy -> ref x, read y (one base each)
                i += 2;
                out.push(Site {
                    rec_q: q,
                    ref_pos: r,
                });
                q += 1;
                r += 1;
            }
            b'=' => {
                // =SEQ -> matched bases
                let start = i;
                while i < cs.len() && cs[i].is_ascii_alphabetic() {
                    i += 1;
                }
                let n = i - start;
                q += n;
                r += n;
            }
            b'+' => {
                // +seq -> insertion (query only)
                let start = i;
                while i < cs.len() && cs[i].is_ascii_alphabetic() {
                    i += 1;
                }
                q += i - start;
            }
            b'-' => {
                // -seq -> deletion (ref only)
                let start = i;
                while i < cs.len() && cs[i].is_ascii_alphabetic() {
                    i += 1;
                }
                r += i - start;
            }
            b'~' => {
                // ~ss<int>ss -> intron: 2 donor bases, int, 2 acceptor bases; ref skip.
                i += 2; // donor
                let mut n = 0usize;
                while i < cs.len() && cs[i].is_ascii_digit() {
                    n = n * 10 + (cs[i] - b'0') as usize;
                    i += 1;
                }
                i += 2; // acceptor
                r += n;
            }
            _ => break, // unknown token; stop parsing defensively
        }
    }
    out
}

/// Parse an `MD:Z` string together with the CIGAR. MD describes matches,
/// substitutions, and deletions along the M/=/X columns of the alignment.
///
/// MD tokens: a run of digits = that many matched columns, a letter = one
/// substituted reference base (one column), `^SEQ` = a deletion (ref only,
/// aligns with a CIGAR D op). Insertions and clips do not appear in MD, so a
/// single MD match-run may span multiple CIGAR M ops separated by I — hence the
/// match counter (`run`) persists across ops.
pub fn from_md(md: &[u8], ops: &CigarOps, ref_start: usize) -> Vec<Site> {
    let mut out = Vec::new();
    let mut q = 0usize;
    let mut r = ref_start;
    let mut mi = 0usize; // index into md
    let mut run = 0usize; // remaining matched columns from the current digit token

    // Consume one aligned column: returns true if it is a substitution.
    let next_column = |md: &[u8], mi: &mut usize, run: &mut usize| -> Option<bool> {
        loop {
            if *run > 0 {
                *run -= 1;
                return Some(false); // match column
            }
            if *mi >= md.len() {
                return None;
            }
            let b = md[*mi];
            if b.is_ascii_digit() {
                let mut n = 0usize;
                while *mi < md.len() && md[*mi].is_ascii_digit() {
                    n = n * 10 + (md[*mi] - b'0') as usize;
                    *mi += 1;
                }
                *run = n;
                continue;
            } else if b == b'^' {
                // deletion token: not a query/M column, skip it here.
                *mi += 1;
                while *mi < md.len() && md[*mi].is_ascii_alphabetic() {
                    *mi += 1;
                }
                continue;
            } else if b.is_ascii_alphabetic() {
                *mi += 1;
                return Some(true); // substitution column
            } else {
                *mi += 1; // stray separator
            }
        }
    };

    for &(len, k) in ops {
        match k {
            b'M' | b'=' | b'X' => {
                for i in 0..len {
                    match next_column(md, &mut mi, &mut run) {
                        Some(true) => out.push(Site {
                            rec_q: q + i,
                            ref_pos: r + i,
                        }),
                        Some(false) => {}
                        None => break,
                    }
                }
                q += len;
                r += len;
            }
            b'I' | b'S' => q += len,
            b'D' | b'N' => r += len,
            _ => {}
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use crate::align::parse_cigar_str;

    fn qpos(sites: &[super::Site]) -> Vec<usize> {
        sites.iter().map(|s| s.rec_q).collect()
    }

    #[test]
    fn cs_short_form_substitutions() {
        // lead_soft=5, ref_start=100. cs: 3 match, sub, 2 match, sub.
        let sites = super::from_cs(b":3*ac:2*gt", 5, 100);
        // subs at SEQ positions 5+3=8 and 8+1+2=11.
        assert_eq!(qpos(&sites), vec![8, 11]);
        assert_eq!(sites[0].ref_pos, 103);
        assert_eq!(sites[1].ref_pos, 106);
    }

    #[test]
    fn cs_handles_indels() {
        // insertion (+ query only) and deletion (- ref only) shift coords.
        let sites = super::from_cs(b":2+ac:2*gt", 0, 0);
        // after :2 -> q=2,r=2; +ac -> q=4; :2 -> q=6,r=4; *gt at q=6,r=4
        assert_eq!(qpos(&sites), vec![6]);
        assert_eq!(sites[0].ref_pos, 4);
    }

    #[test]
    fn md_with_insertion_spanning_run() {
        // CIGAR 5M2I5M, MD "3A6": one substitution 4 bases into the alignment,
        // where the MD match-run (6) spans the I between the two M ops.
        let ops = parse_cigar_str(b"5M2I5M").unwrap();
        let sites = super::from_md(b"3A6", &ops, 0);
        // 3 matches, sub at column 3 (rec_q 3), then 6 matches across the 2M ops.
        assert_eq!(qpos(&sites), vec![3]);
        assert_eq!(sites[0].ref_pos, 3);
    }

    #[test]
    fn cigar_x_mismatches() {
        let ops = parse_cigar_str(b"3M1X4M").unwrap();
        let sites = super::from_cigar_x(&ops);
        assert_eq!(qpos(&sites), vec![3]);
    }

    #[test]
    fn reference_comparison() {
        let ops = parse_cigar_str(b"5M").unwrap();
        let seq = b"ACGTA";
        let refseq = b"ACCTA"; // differ at index 2
        let sites = super::from_reference(&ops, seq, 0, refseq);
        assert_eq!(qpos(&sites), vec![2]);
        assert_eq!(sites[0].ref_pos, 2);
    }
}
