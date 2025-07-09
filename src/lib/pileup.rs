use crate::{modtype::ModType, strand::Strand};
use ahash::AHashMap as HashMap;
use anyhow::anyhow;
use anyhow::Result;
use strum::IntoEnumIterator;
use strum_macros::{EnumIter, ToString};
use atoi;
use csv::{ByteRecord, ReaderBuilder};
use log::{debug, info, warn};
use serde::Serialize;
use std::collections::VecDeque;
use std::io::Read as IoRead;
use std::str::FromStr;
use rust_htslib::tbx::{self, Read};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, EnumIter)]
pub enum MethylBedFields {
    Chrom,
    Start,
    End,
    ModCode,
    Score,
    Strand,
    StartCompat,
    EndCompat,
    Color,
    NValidCov,
    FractionModified,
    NMod,
    NCanonical,
    NOtherMod,
    NDelete,
    NFail,
    NDiff,
    NNocall, // Placeholder for future use
}

impl MethylBedFields {
    pub fn all_as_string() -> Vec<String> {
        MethylBedFields::iter().map(|f| f.to_string()).collect()
    }
}

impl FromStr for MethylBedFields {
    type Err = anyhow::Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "chrom" => Ok(MethylBedFields::Chrom),
            "start_position" => Ok(MethylBedFields::Start),
            "end_position" => Ok(MethylBedFields::End),
            "modified_base_code" => Ok(MethylBedFields::ModCode),
            "score" => Ok(MethylBedFields::Score),
            "strand" => Ok(MethylBedFields::Strand),
            "start_position_compat" => Ok(MethylBedFields::StartCompat),
            "end_position_compat" => Ok(MethylBedFields::EndCompat),
            "color" => Ok(MethylBedFields::Color),
            "Nvalid_cov" => Ok(MethylBedFields::NValidCov),
            "fraction_modified" => Ok(MethylBedFields::FractionModified),
            "Nmod" => Ok(MethylBedFields::NMod),
            "Ncanonical" => Ok(MethylBedFields::NCanonical),
            "Nother_mod" => Ok(MethylBedFields::NOtherMod),
            "Ndelete" => Ok(MethylBedFields::NDelete),
            "Nfail" => Ok(MethylBedFields::NFail),
            "Ndiff" => Ok(MethylBedFields::NDiff),
            "Nnocall" => Ok(MethylBedFields::NNocall), // Placeholder for future use
            _ => Err(anyhow!("Unknown field: {}", s)),
        }
    }
}

impl ToString for MethylBedFields {
    fn to_string(&self) -> String {
        match self {
            MethylBedFields::Chrom => "chrom".to_string(),
            MethylBedFields::Start => "start_position".to_string(),
            MethylBedFields::End => "end_position".to_string(),
            MethylBedFields::ModCode => "modified_base_code".to_string(),
            MethylBedFields::Score => "score".to_string(),
            MethylBedFields::Strand => "strand".to_string(),
            MethylBedFields::StartCompat => "start_position_compat".to_string(),
            MethylBedFields::EndCompat => "end_position_compat".to_string(),
            MethylBedFields::Color => "color".to_string(),
            MethylBedFields::NValidCov => "Nvalid_cov".to_string(),
            MethylBedFields::FractionModified => "fraction_modified".to_string(),
            MethylBedFields::NMod => "Nmod".to_string(),
            MethylBedFields::NCanonical => "Ncanonical".to_string(),
            MethylBedFields::NOtherMod => "Nother_mod".to_string(),
            MethylBedFields::NDelete => "Ndelete".to_string(),
            MethylBedFields::NFail => "Nfail".to_string(),
            MethylBedFields::NDiff => "Ndiff".to_string(),
            MethylBedFields::NNocall => "Nnocall".to_string(), // Placeholder for future use

        }
    }
}

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

#[derive(Debug, Clone, PartialEq, Eq, Serialize)]
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

impl PileupRecord {
    pub fn to_bed_fields(&self) -> Result<Vec<String>> {
        let chrom = self.reference.clone();
        let start = self.position;
        let end = start + 1;
        let mod_code = self.mod_type.to_pileup_code().to_owned();
        let score = self.n_valid_cov;
        let strand = self.strand.to_string();

        let fraction_modified = if self.n_valid_cov > 0 {
            format!("{:.1}", 100.0 * self.n_mod as f32 / self.n_valid_cov as f32)
        } else {
            "0.0".to_string()
        };
        let bed_fields = vec![
            chrom,
            start.to_string(),
            end.to_string(),
            mod_code,
            score.to_string(),
            strand,
            start.to_string(),     // compatibility
            end.to_string(),       // compatibility
            "255,0,0".to_string(), // BED RGB
            self.n_valid_cov.to_string(),
            fraction_modified,
            self.n_mod.to_string(),
            self.n_canonical.to_string(),
            "NA".to_string(),
            "NA".to_string(),
            "NA".to_string(),
            self.n_diff.to_string(),
            "NA".to_string(),
        ];
        Ok(bed_fields)
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PileupChunk {
    pub reference: String,
    pub records: Vec<PileupRecord>,
}

impl PileupChunk {
    pub fn new(reference: String, records: Vec<PileupRecord>) -> Self {
        Self {
            reference,
            records,
        }
    }

    pub fn is_empty(&self) -> bool {
        self.records.is_empty()
    }


    pub fn filter(&self, filter_fn: impl Fn(&PileupRecord) -> bool) -> Self {
        let filtered_records: Vec<PileupRecord> = self
            .records
            .iter()
            .filter(|r| filter_fn(r))
            .cloned()
            .collect();
        Self {
            reference: self.reference.clone(),
            records: filtered_records,
        }
    }
    pub fn filter_min_cov(
        &self,
        min_cov: u32,
    ) -> Self {
        self.filter(|r| r.n_valid_cov >= min_cov)
    }

    pub fn filter_fraction_n_diff(
        &self,
        fraction: f32,
    ) -> Self {
        self.filter(|r| {
            if r.n_valid_cov == 0 {
                return false; // Avoid division by zero
            }
            let fraction_n_diff = r.n_diff as f32 / (r.n_diff + r.n_valid_cov) as f32;
            fraction_n_diff >= fraction
        })
    }
}

pub struct PileupChunkReader<R: IoRead> {
    reader: csv::Reader<R>,
    buffer: VecDeque<ByteRecord>,
    min_cov: u32,
    pub eof_reached: bool,
}

impl<R: IoRead> PileupChunkReader<R> {
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

            if let Ok(parsed_record) =
                parse_and_validate_pileup_record(&record, self.min_cov, &field_mapping)
            {
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

            if let Ok(parsed_record) =
                parse_and_validate_pileup_record(&record, self.min_cov, &field_mapping)
            {
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



pub struct PileupTabixReader {
    pub reader: tbx::Reader,
    pub min_cov: u32,
    pub threads: usize,
}

impl PileupTabixReader {
    pub fn new(inner: String, min_cov: u32, threads: usize) -> Result<Self> {
        let mut reader = tbx::Reader::from_path(inner)?;
        reader.set_threads(threads)
            .map_err(|e| anyhow!("Failed to set threads for tbx reader: {}", e))?;
        Ok(Self {
            reader,
            min_cov,
            threads,
        })
    }
    pub fn fetch_contig(&mut self, contig: &str) -> Result<(), Box<dyn std::error::Error>> {
        let tid = match self.reader.tid(contig) {
            Ok(tid) => tid,
            Err(_) => panic!("Could not resolve contig name '{}'", contig),
        };
        self.reader.fetch(tid, 0, i32::MAX as u64)
            .expect("Could not seek");
        Ok(())
    }

    pub fn read_all_from_contig(&mut self) -> Option<PileupChunk> {
        let mut parsed_records = Vec::with_capacity(1_000);
        let mut buf = Vec::new();
        
        let field_mapping = FieldMapping::default();
        while let Ok(line) = self.reader.read(&mut buf) {
            if !line {
                break;
            } 
            let fields: Vec<&[u8]> = buf.split(|b| *b == b'\t').collect(); 
            match parse_and_validate_pileup_record_bytes(&fields, self.min_cov, &field_mapping) {
                Ok(parsed_record) => parsed_records.push(parsed_record),
                Err(e) => continue, // Skip invalid records
            }
        }

        if parsed_records.is_empty() {
            return None;
        }

        Some(PileupChunk {
            reference: parsed_records[0].reference.clone(),
            records: parsed_records,
        })
    }
}



pub fn parse_and_validate_pileup_record_bytes(
    fields: &[&[u8]],
    min_cov: u32,
    mapping: &FieldMapping,
) -> Result<PileupRecord> {
    let nvalid_idx = mapping.idx(PileupField::NValidCov).ok_or_else(|| anyhow!("Missing NValidCov field"))?;
    let n_valid_cov = atoi::atoi(fields[nvalid_idx]).ok_or_else(|| anyhow!("Invalid n_valid_cov"))?;
    if n_valid_cov < min_cov {
        return Err(anyhow!("Coverage below minimum threshold"));
    }
    let ref_idx = mapping.idx(PileupField::Reference).ok_or_else(|| anyhow!("Missing Reference field"))?;
    let pos_idx = mapping.idx(PileupField::Position).ok_or_else(|| anyhow!("Missing Position field"))?;
    let strand_idx = mapping.idx(PileupField::Strand).ok_or_else(|| anyhow!("Missing Strand field"))?;
    let modtype_idx = mapping.idx(PileupField::ModType).ok_or_else(|| anyhow!("Missing ModType field"))?;
    let nmod_idx = mapping.idx(PileupField::NMod).ok_or_else(|| anyhow!("Missing NMod field"))?;
    let ncan_idx = mapping.idx(PileupField::NCanonical).ok_or_else(|| anyhow!("Missing NCanonical field"))?;
    let ndiff_idx = mapping.idx(PileupField::NDiff).ok_or_else(|| anyhow!("Missing NDiff field"))?;

    let reference = std::str::from_utf8(fields[ref_idx])?.to_string();
    let position = atoi::atoi(fields[pos_idx]).ok_or_else(|| anyhow!("Invalid position"))?;
    let strand_str = std::str::from_utf8(fields[strand_idx])?;
    let strand = strand_str.parse::<Strand>()
        .map_err(|_| anyhow!("Could not parse pileup strand value"))?;
    let modtype_str = std::str::from_utf8(fields[modtype_idx])?;
    let mod_type = modtype_str.parse::<ModType>()
        .map_err(|_| anyhow!("Could not parse pileup mod_type value"))?;



    let n_mod = atoi::atoi(fields[nmod_idx]).ok_or_else(|| anyhow!("Invalid n_mod"))?;
    let n_canonical = atoi::atoi(fields[ncan_idx]).ok_or_else(|| anyhow!("Invalid n_canonical"))?;
    let n_diff = atoi::atoi(fields[ndiff_idx]).ok_or_else(|| anyhow!("Invalid n_diff"))?;

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

fn line_to_byterecord(line: &[u8]) -> ByteRecord {
    let mut record = ByteRecord::new();
    record.extend(line.split(|b| *b == b'\t'));
    record
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
    fn test_pileup_record_to_bed_fields() {
        let record = PileupRecord {
            reference: "contig_1".to_string(),
            position: 10,
            strand: Strand::Positive,
            mod_type: ModType::SixMA,
            n_mod: 1,
            n_valid_cov: 10,
            n_canonical: 9,
            n_diff: 3,
        };
        let bed_fields = record.to_bed_fields().unwrap();
        assert_eq!(bed_fields.len(), 18);
        assert_eq!(bed_fields[0], "contig_1");
        assert_eq!(bed_fields[1], "10");    
        assert_eq!(bed_fields[2], "11"); // end position is start + 1
        assert_eq!(bed_fields[3], "a");
        assert_eq!(bed_fields[4], "10");
        assert_eq!(bed_fields[5], "+");
        assert_eq!(bed_fields[6], "10"); // start position compatibility
        assert_eq!(bed_fields[7], "11"); // end position compatibility
        assert_eq!(bed_fields[8], "255,0,0"); // BED RGB
        assert_eq!(bed_fields[9], "10"); // Nvalid_cov
        assert_eq!(bed_fields[10], "10.0"); // fraction_modified
        assert_eq!(bed_fields[11], "1"); // Nmod
        assert_eq!(bed_fields[12], "9"); // Ncanonical
        assert_eq!(bed_fields[13], "NA"); // Nother_mod
        assert_eq!(bed_fields[14], "NA"); // Ndelete
        assert_eq!(bed_fields[15], "NA"); // Nfail
        assert_eq!(bed_fields[16], "3"); // Ndiff
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
