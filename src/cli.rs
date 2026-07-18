use clap::Parser;
use std::path::PathBuf;

/// Extract bad read alignments (long clips and mismatch-dense regions) from a BAM.
#[derive(Parser, Debug)]
#[command(name = "badalign", version, about)]
pub struct Args {
    /// Input BAM file.
    pub bam: PathBuf,

    /// Reference genome FASTA (loaded fully into RAM).
    #[arg(short = 'r', long = "ref", value_name = "FILE")]
    pub reference: Option<PathBuf>,

    /// Length of flanking regions around a clip / dense region.
    #[arg(short = 'l', value_name = "INT", default_value_t = 250)]
    pub flank: usize,

    /// Minimum base quality for a mismatch to count (type-2).
    #[arg(short = 'Q', value_name = "INT", default_value_t = 20)]
    pub min_baseq: u8,

    /// Window length to scan for a dense region (type-2).
    #[arg(short = 'w', value_name = "INT", default_value_t = 100)]
    pub window: usize,

    /// Minimum clip length to report (type-1).
    #[arg(short = 'c', value_name = "INT", default_value_t = 40)]
    pub min_clip: usize,

    /// Minimum number of high-quality mismatches in a window (type-2).
    #[arg(short = 'm', value_name = "INT", default_value_t = 10)]
    pub min_mismatch: usize,
}
