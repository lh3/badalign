//! Whole-genome reference held in RAM (no faidx), keyed by contig name.

use anyhow::{Context, Result};
use std::collections::HashMap;
use std::path::Path;

pub struct RefGenome {
    seqs: HashMap<String, Vec<u8>>,
}

impl RefGenome {
    pub fn load(path: &Path) -> Result<Self> {
        let mut reader = noodles::fasta::io::reader::Builder
            .build_from_path(path)
            .with_context(|| format!("opening reference {}", path.display()))?;
        let mut seqs = HashMap::new();
        for result in reader.records() {
            let record = result.context("reading FASTA record")?;
            let name = String::from_utf8_lossy(record.name()).into_owned();
            seqs.insert(name, record.sequence().as_ref().to_vec());
        }
        Ok(Self { seqs })
    }

    pub fn get(&self, ctg: &str) -> Option<&[u8]> {
        self.seqs.get(ctg).map(|v| v.as_slice())
    }
}
