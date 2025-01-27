use anyhow::anyhow;
use anyhow::Result;
use atoi;
use csv::{ByteRecord, ReaderBuilder};
use log::{debug, info, warn};
use crate::{
    modtype::ModType, 
    strand::Strand
};
use std::collections::VecDeque;
use std::io::Read;
use ahash::AHashMap as HashMap;


#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum PileupField {
    Reference,
    Position,
    Strand,
    ModType,
    NMod,
    NValidCov,
    NCanonical,
    NDiff,
}

#[derive(Debug, Clone)]
pub struct FieldMapping {
    pub mapping: HashMap<PileupField, usize>,
}

impl FieldMapping {
    pub fn new() -> Self {
        Self {
            mapping: HashMap::new(),
        }
    }

    pub fn with_field(mut self, field: PileupField, idx: usize) -> Self {
        self.mapping.insert(field, idx);
        self
    }

    pub fn idx(&self, field: PileupField) -> Option<usize> {
        self.mapping.get(&field).copied()
    }
}

impl Default for FieldMapping {
    fn default() -> Self {
        Self::new()
        .with_field(PileupField::Reference, 0)
        .with_field(PileupField::Position, 1)
        .with_field(PileupField::Strand, 5)
        .with_field(PileupField::ModType, 3)
        .with_field(PileupField::NMod, 11)
        .with_field(PileupField::NValidCov, 9)
        .with_field(PileupField::NCanonical, 12)
        .with_field(PileupField::NDiff, 16)
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PileupRecord {
    pub reference: String,
    pub position: usize,
    pub strand: Strand,
    pub mod_type: ModType,
    pub n_mod: u32,
    pub n_valid_cov: u32,
    pub n_canonical: u32,
    pub n_diff: u32,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PileupChunk {
    pub reference: String,
    pub records: Vec<PileupRecord>,
}

pub struct PileupChunkReader<R: Read> {
    reader: csv::Reader<R>,
    buffer: VecDeque<ByteRecord>,
    min_cov: u32,
    pub eof_reached: bool,
}

impl<R: Read> PileupChunkReader<R> {
    pub fn new(inner: R, min_cov: u32) -> Self {
        let reader = ReaderBuilder::new()
            .delimiter(b'\t')
            .has_headers(false)
            .buffer_capacity(128 * (1 << 10))
            .from_reader(inner);
        Self {
            reader,
            buffer: VecDeque::new(),
            min_cov,
            eof_reached: false,
        }
    }

    /// Reads the next chunk of records grouped by the same reference
    pub fn next_chunk(&mut self) -> Option<PileupChunk> {
        let mut parsed_records = Vec::new();
        let mut current_reference = None;
        let mut record = ByteRecord::new();
        let field_mapping = FieldMapping::default();
        // Process records from the buffer
        while let Some(record) = self.buffer.pop_front() {
            let reference = std::str::from_utf8(record.get(0).unwrap_or(b""))
                .unwrap_or("")
                .to_string();

            match &current_reference {
                Some(current_ref) if reference != *current_ref => {
                    self.buffer.push_front(record);
                    break;
                }
                None => {
                    current_reference = Some(reference.clone());
                }
                _ => {}
            }

            if let Ok(parsed_record) = parse_and_validate_pileup_record(&record, self.min_cov, &field_mapping) {
                parsed_records.push(parsed_record);
            }
        }

        // Load the next batch of records
        while let Ok(has_record) = self.reader.read_byte_record(&mut record) {
            if !has_record {
                self.eof_reached = true;
                break;
            }
            let n_valid_cov: u32 = atoi::atoi(record.get(9).expect("Failed to get n_valid_cov"))
                .expect("Failed to parse n_valid_cov value");
            if n_valid_cov < self.min_cov {
                continue;
            }
            let reference = std::str::from_utf8(record.get(0).unwrap_or(b"")).unwrap();

            match current_reference {
                Some(current_ref) if reference != current_ref => {
                    self.buffer.push_front(record);
                    break;
                }
                None => {
                    current_reference = Some(reference.to_string());
                }
                _ => {}
            }

            if let Ok(parsed_record) = parse_and_validate_pileup_record(&record, self.min_cov, &field_mapping) {
                parsed_records.push(parsed_record);
            }
        }
        if parsed_records.is_empty() {
            return None;
        }
        // Get the reference from the first record
        let reference_out = parsed_records[0].reference.clone();
        Some(PileupChunk {
            reference: reference_out,
            records: parsed_records,
        })
    }

    pub fn load_n_chunks(&mut self, n: usize) -> Option<Vec<PileupChunk>> {
        let mut chunks = Vec::new();
        for _ in 0..n {
            if let Some(chunk) = self.next_chunk() {
                chunks.push(chunk);
            } else {
                break; // Stop if no more chunks are available
            }
        }

        if chunks.is_empty() {
            None
        } else {
            Some(chunks)
        }
    }
}
fn parse_and_validate_pileup_record(
    record: &ByteRecord,
    min_cov: u32,
    field_mapping: &FieldMapping,
) -> Result<PileupRecord> {
    let reference_idx = field_mapping
        .idx(PileupField::Reference)
        .ok_or_else(|| anyhow!("No column mapping for Reference"))?;
    let position_idx = field_mapping
        .idx(PileupField::Position)
        .ok_or_else(|| anyhow!("No column mapping for Position"))?;
    let mod_type_idx = field_mapping
        .idx(PileupField::ModType)
        .ok_or_else(|| anyhow!("No column mapping for ModType"))?;
    let strand_idx = field_mapping
        .idx(PileupField::Strand)
        .ok_or_else(|| anyhow!("No column mapping for Strand"))?;
    let n_valid_cov_idx = field_mapping
        .idx(PileupField::NValidCov)
        .ok_or_else(|| anyhow!("No column mapping for NValidCov"))?;
    let n_mod_idx = field_mapping
        .idx(PileupField::NMod)
        .ok_or_else(|| anyhow!("No column mapping for NMod"))?;
    let n_canonical_idx = field_mapping
        .idx(PileupField::NCanonical)
        .ok_or_else(|| anyhow!("No column mapping for NCanonical"))?;
    let n_diff_idx = field_mapping
        .idx(PileupField::NDiff)
        .ok_or_else(|| anyhow!("No column mapping for NDiff"))?;

    let reference = std::str::from_utf8(record.get(reference_idx).unwrap_or(b""))?.to_string();
    let position = atoi::atoi::<usize>(record.get(position_idx).unwrap_or(b""))
        .ok_or_else(|| anyhow!("Could not parse pileup position"))?;

    let strand = std::str::from_utf8(record.get(strand_idx).unwrap_or(b""))?
        .parse::<Strand>()
        .map_err(|_| anyhow!("Could not parse pileup strand value"))?;

    let mod_type = std::str::from_utf8(record.get(mod_type_idx).unwrap_or(b""))?
        .parse::<ModType>()
        .map_err(|_| anyhow!("Could not parse pileup mod_type value"))?;

    let n_valid_cov = atoi::atoi::<u32>(record.get(n_valid_cov_idx).unwrap_or(b""))
        .ok_or_else(|| anyhow!("Invalid n_valid_cov value"))?;

    // coverage check
    if n_valid_cov < min_cov {
        return Err(anyhow!("Coverage below minimum threshold"));
    }

    let n_mod = atoi::atoi::<u32>(record.get(n_mod_idx).unwrap_or(b""))
        .ok_or_else(|| anyhow!("Could not parse pileup n_mod value"))?;
    let n_canonical = atoi::atoi::<u32>(record.get(n_canonical_idx).unwrap_or(b""))
        .ok_or_else(|| anyhow!("Could not parse pileup n_canonical value"))?;
    let n_diff = atoi::atoi::<u32>(record.get(n_diff_idx).unwrap_or(b""))
        .ok_or_else(|| anyhow!("Could not parse pileup n_diff value"))?;

    Ok(PileupRecord {
        reference,
        position,
        strand,
        mod_type,
        n_mod,
        n_valid_cov,
        n_canonical,
        n_diff,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::{Seek, SeekFrom, Write};
    use tempfile::NamedTempFile;

    /// Creates one line of a pileup-like record, with many columns separated by tabs.
    ///
    /// Columns:
    ///   1) reference
    ///   2) position
    ///   3) '.'
    ///   4) mod_type
    ///   5) '.'
    ///   6) strand
    ///   7) '.'
    ///   8) '.'
    ///   9) '.'
    ///   10) n_valid_cov
    ///   11) '.'
    ///   12) n_mod
    ///   13) n_canonical
    ///   14) '.'
    ///   15) '.'
    ///   16) '.'
    ///   17) n_diff
    ///   18) '.'
    ///
    /// Note: Adjust the placeholders or column order to match your actual parsing logic.
    fn create_pileup_line(
        reference: &str,
        position: usize,
        strand: &str,
        mod_type: &str,
        n_mod: u32,
        n_valid_cov: u32,
        n_canonical: u32,
        n_diff: u32,
    ) -> String {
        format!(
            "{}\t{}\t.\t{}\t.\t{}\t.\t.\t.\t{}\t.\t{}\t{}\t.\t.\t.\t{}\t.\n",
            reference,   // 1
            position,    // 2
            mod_type,    // 4
            strand,      // 6
            n_valid_cov, // 10
            n_mod,       // 12
            n_canonical, // 13
            n_diff       // 17
        )
    }
    fn create_temp_file(data: &[u8]) -> NamedTempFile {
        let mut temp_file = NamedTempFile::new().expect("Failed to create temp file");
        temp_file.write_all(data).expect("Failed to write data");
        temp_file
            .seek(SeekFrom::Start(0))
            .expect("Failed to seek to start");
        temp_file
    }

    #[test]
    fn test_new_pileup_record() {
        let record = PileupRecord {
            reference: "contig_1".to_string(),
            position: 10,
            strand: Strand::Positive,
            mod_type: ModType::SixMA,
            n_mod: 1,
            n_valid_cov: 10,
            n_canonical: 2,
            n_diff: 3,
        };
        assert_eq!(record.reference, "contig_1");
        assert_eq!(record.position, 10);
        assert_eq!(record.strand, Strand::Positive);
        assert_eq!(record.mod_type, ModType::SixMA);
        assert_eq!(record.n_mod, 1);
        assert_eq!(record.n_valid_cov, 10);
        assert_eq!(record.n_canonical, 2);
        assert_eq!(record.n_diff, 3);
    }

    #[test]
    fn test_new_pileup_chunk() {
        let record = PileupRecord {
            reference: "contig_1".to_string(),
            position: 10,
            strand: Strand::Positive,
            mod_type: ModType::SixMA,
            n_mod: 1,
            n_valid_cov: 10,
            n_canonical: 2,
            n_diff: 3,
        };
        let chunk = PileupChunk {
            reference: "contig_1".to_string(),
            records: vec![record],
        };
        assert_eq!(chunk.reference, "contig_1");
        assert_eq!(chunk.records.len(), 1);
    }

    #[test]
    fn test_new_pileup_chunk_reader() {
        let data = b"contig_1\t10\t.\ta\t.\t+\t.\t.\t.\t10\t.\t1\t2\t.\t.\t.\t3\t.\n";
        let tempfile = create_temp_file(data);
        let reader = PileupChunkReader::new(tempfile.reopen().unwrap(), 1);
        assert_eq!(reader.min_cov, 1);
    }
    /// Use tempfile to write test data to a file, then test with PileupChunkReader
    #[test]
    fn test_with_tempfile() {
        let mut data = String::new();
        data.push_str(create_pileup_line("contig_1", 0, "+", "a", 10, 10, 0, 0).as_str());
        data.push_str(create_pileup_line("contig_1", 1, "+", "a", 10, 10, 0, 0).as_str());
        data.push_str(create_pileup_line("contig_2", 0, "+", "m", 10, 10, 0, 0).as_str());

        let tempfile = create_temp_file(data.as_bytes());
        let mut reader = PileupChunkReader::new(tempfile.reopen().unwrap(), 1);

        // First chunk should correspond to 'contig_1'
        let chunk_opt = reader.next_chunk();
        assert!(chunk_opt.is_some());
        let chunk = chunk_opt.unwrap();
        assert_eq!(chunk.reference, "contig_1");
        assert_eq!(chunk.records.len(), 2);

        // Next chunk should correspond to 'contig_2'
        let chunk_opt_2 = reader.next_chunk();
        assert!(chunk_opt_2.is_some());
        let chunk_2 = chunk_opt_2.unwrap();
        assert_eq!(chunk_2.reference, "contig_2");
        // Only one record in the sample data for contig_2
        assert_eq!(chunk_2.records.len(), 1);

        // No more data, so the next should be None
        let chunk_opt_3 = reader.next_chunk();
        assert!(chunk_opt_3.is_none());
        assert!(reader.eof_reached);
    }

    #[test]
    fn test_low_coverage_filtering() {
        let mut data = String::new();
        data.push_str(create_pileup_line("contig_1", 0, "+", "a", 81, 0, 0, 0).as_str());
        data.push_str(create_pileup_line("contig_1", 1, "+", "a", 82, 3, 0, 0).as_str());
        data.push_str(create_pileup_line("contig_2", 0, "+", "m", 84, 10, 0, 0).as_str());

        let tempfile = create_temp_file(data.as_bytes());
        let mut reader = PileupChunkReader::new(tempfile.reopen().unwrap(), 2);

        // First chunk should correspond to 'contig_2'
        let chunk_opt = reader.next_chunk();
        assert!(chunk_opt.is_some());
        let chunk = chunk_opt.unwrap();
        assert_eq!(chunk.reference, "contig_1");
        assert_eq!(chunk.records.len(), 1);

        // Next chunk should correspond to 'contig_2'
        let chunk_opt_2 = reader.next_chunk();
        assert!(chunk_opt_2.is_some());
        let chunk_2 = chunk_opt_2.unwrap();
        assert_eq!(chunk_2.reference, "contig_2");
        // Only one record in the sample data for contig_2
        assert_eq!(chunk_2.records.len(), 1);
    }

    #[test]
    fn test_n_chunks() {
        let mut data = String::new();
        data.push_str(create_pileup_line("contig_1", 0, "+", "a", 10, 10, 0, 0).as_str());
        data.push_str(create_pileup_line("contig_2", 1, "+", "a", 10, 10, 0, 0).as_str());
        data.push_str(create_pileup_line("contig_3", 0, "+", "m", 10, 10, 0, 0).as_str());

        let tempfile = create_temp_file(data.as_bytes());
        let mut reader = PileupChunkReader::new(tempfile.reopen().unwrap(), 1);

        let chunks = reader.load_n_chunks(2).unwrap();
        assert_eq!(chunks.len(), 2);
        assert_eq!(chunks[0].reference, "contig_1");
        assert_eq!(chunks[0].records.len(), 1);
        assert_eq!(chunks[1].reference, "contig_2");
        assert_eq!(chunks[1].records.len(), 1);
    }
    #[test]
    fn test_parse_and_validate_record() {
        let record = ByteRecord::from(vec![
            b"contig_1".to_vec(),
            b"6".to_vec(),
            b".".to_vec(),
            b"a".to_vec(),
            b".".to_vec(),
            b"+".to_vec(),
            b".".to_vec(),
            b".".to_vec(),
            b".".to_vec(),
            b"3".to_vec(),
            b".".to_vec(),
            b"1".to_vec(),
            b"2".to_vec(),
            b".".to_vec(),
            b".".to_vec(),
            b".".to_vec(),
            b"4".to_vec(),
            b".".to_vec(),
        ]);
        let field_mapping = FieldMapping::default();
        let parsed_record = parse_and_validate_pileup_record(&record, 1, &field_mapping).unwrap();
        assert_eq!(parsed_record.reference, "contig_1");
        assert_eq!(parsed_record.position, 6);
        assert_eq!(parsed_record.strand, Strand::Positive);
        assert_eq!(parsed_record.mod_type, ModType::SixMA);
        assert_eq!(parsed_record.n_mod, 1);
        assert_eq!(parsed_record.n_valid_cov, 3);
        assert_eq!(parsed_record.n_canonical, 2);
        assert_eq!(parsed_record.n_diff, 4);
    }
}
