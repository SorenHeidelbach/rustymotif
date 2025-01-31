use crate::sequence::Contig;
use crate::sequence::MethylationThresholds;
use anyhow::Context;
use anyhow::Result;
use seq_io::fasta::{Reader, Record};
use std::collections::HashMap;
use std::path::Path;

pub fn read_fasta_file(file_path: &Path) -> Result<HashMap<String, String>> {
    let mut reader = Reader::from_path(file_path)
        .with_context(|| format!("Error reading fasta file: {}", file_path.display()))?;

    let mut records = HashMap::new();

    while let Some(record_result) = reader.next() {
        let record = record_result.with_context(|| "Error reading fasta record")?;
        let id = record
            .id()
            .map(|id| id.to_string())
            .with_context(|| "Error getting fasta record id")?;
        let sequence = String::from_utf8(record.owned_seq())
            .with_context(|| "Error getting fasta record sequence")?;

        records.insert(id, sequence);
    }
    Ok(records)
}
