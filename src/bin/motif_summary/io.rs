use std::collections::HashMap;
use std::path::Path;
use std::fs::File;
use std::fs;
use std::io::{BufReader, BufRead};
use anyhow::Result;
use csv::{WriterBuilder, Writer};
use rustymotif_utils::motif::Motif;
use log::{info, debug, warn};

use crate::data::AnnotatedPileupRecord;
use crate::data;

pub fn load_contig_bin(
    contig_bin_path: &Path,
) -> Result<HashMap<String, Vec<String>>> {
    let file = File::open(contig_bin_path)
        .expect("Error opening contig bin file");
    
    let reader = BufReader::new(file);
    let mut bin_contigs = HashMap::new();
    for line in reader.lines() {
        let line = line.expect("Error reading line from contig bin file");
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 2 {
            continue; // Skip lines that do not have at least two parts
        }
        let contig = parts[0].to_string();
        let bin  = parts[1].to_string();
        bin_contigs
            .entry(bin)
            .or_insert_with(Vec::new)
            .push(contig);
    }
    Ok(bin_contigs)


}
pub fn load_bin_motifs(
    bin_motifs_path: &Path,
) -> Result<HashMap<String, Vec<Motif>>> {
    let file = File::open(bin_motifs_path)?;
    
    let reader = BufReader::new(file);
    let mut bin_motifs = HashMap::new();
    for line in reader.lines().skip(1) {
        let line = line?;
        let parts: Vec<&str> = line.split('\t').collect();
        let bin = parts[0].to_string();


        // Forward motif
        let motif_seq = parts[1];
        let motif_position = parts[2].parse::<u8>()?;
        let motif_modtype = parts[3];
        let motif = Motif::new(motif_seq, motif_modtype, motif_position)?;
        bin_motifs
                .entry(bin.clone())
                .or_insert_with(Vec::new)
                .push(motif.clone());
        // Reverse motif
        if parts[7].is_empty() {
            // If there is no reverse motif, we can skip this line
            
            continue;
        }

        let reverse_motif_seq = parts[7];
        let reverse_motif_position = parts[8].parse::<u8>()?;
        let reverse_motif = Motif::new(reverse_motif_seq, motif_modtype, reverse_motif_position)?;

        if reverse_motif.eq(&motif) {
            continue;
        }

        bin_motifs
            .entry(bin)
            .or_insert_with(Vec::new)
            .push(reverse_motif);
    }
    // Ensure all motifs are unique within each bin
    for motifs in bin_motifs.values_mut() {
        let mut unique_motifs = Vec::new();
        let mut seen = HashMap::new();
        for motif in motifs.iter() {
            if !seen.contains_key(&motif.as_pretty_string()) {
                seen.insert(motif.as_pretty_string(), true);
                unique_motifs.push(motif.clone());
            }
        }
        *motifs = unique_motifs;
    }
    Ok(bin_motifs)
}





pub struct AnnotatedRecordWriter {
    file: File,
    writer: Writer<File>,
}

impl AnnotatedRecordWriter {
    pub fn new(file_path: &Path) -> Result<Self> {
        let file = File::create(file_path)?;
        let writer = WriterBuilder::new()
            .delimiter(b'\t')
            .has_headers(true)
            .from_writer(file.try_clone()?);
        Ok(Self { file, writer })
    }

    pub fn write_header(&mut self) -> Result<()> {
        self.writer.write_record(&[
            "contig_id",
            "mod_type",
            "position",
            "strand",
            "n_mod",
            "n_canonical",
            "n_valid_cov",
            "n_diff",                    
            "in_motif",
            "motifs",
        ])?;
        self.writer.flush()?;
        Ok(())
    }

    pub fn write_record(&mut self, record: &AnnotatedPileupRecord) -> Result<()> {
        self.writer.write_record(&[
            &record.pileup_record.reference,
            &record.pileup_record.mod_type.to_string().to_string(),
            &record.pileup_record.position.to_string(),
            &record.pileup_record.strand.to_string(),
            &record.pileup_record.n_mod.to_string(),
            &record.pileup_record.n_canonical.to_string(),
            &record.pileup_record.n_valid_cov.to_string(),
            &record.pileup_record.n_diff.to_string(),
            &record.in_motif.to_string(),
            &record.motifs.iter().map(|m| m.as_pretty_string()).collect::<Vec<String>>().join(","),
        ])?;
        Ok(())
    }

    pub fn flush(&mut self) -> Result<()> {
        self.writer.flush()?;
        self.file.sync_all()?;
        Ok(())
    }

    pub fn write_records_iter<'a, I>(
        &mut self,
        records: I,
    ) -> Result<()> 
    where
        I: IntoIterator<Item = &'a AnnotatedPileupRecord>,
    {
        for record in records {
            self.write_record(record)?;
        }
        self.flush()?;
        Ok(())
    }

}

pub struct AnnotatedRecordSummaryWriter {
    file: File,
    writer: Writer<File>,
} 

impl AnnotatedRecordSummaryWriter {
    pub fn new(path: &Path) -> Result<Self> {
        let file = File::create(path)?;
        let writer = WriterBuilder::new()
            .delimiter(b'\t')
            .has_headers(true)
            .from_writer(file.try_clone()?);
        Ok(Self { file, writer })
    }

    pub fn write_header(&mut self) -> Result<()> {
        self.writer.write_record(&[
            "contig_id",
            "mod_type",
            "strand",
            "motifs",
            "n_motifs",
            "n_records",
            "n_records_below_cov",
            "n_records_in_motif",
            "n_records_not_in_motif",
            "n_mod_in_motif",
            "n_canonical_in_motif",
            "n_valid_cov_in_motif",
            "n_diff_in_motif",
            "n_gt_threshold_in_motif",
            "n_mod_not_in_motif",
            "n_canonical_not_in_motif",
            "n_valid_cov_not_in_motif",
            "n_diff_not_in_motif",
            "n_gt_threshold_not_in_motif",
        ])?;
        self.writer.flush()?;
        Ok(())
    }

    pub fn write_summary(
        &mut self,
        summary: data::AnnotatedPileupSummary,
    ) -> Result<()> {
        for ((mod_type, strand), data) in summary.summary.into_iter() {
            self.writer.write_record(data)?;
        }
        Ok(())
    }

}





#[cfg(test)]
mod tests {
    use super::*;
    use crate::data::AnnotatedPileupRecord;
    use crate::pileup::PileupRecord;
    use rustymotif_utils::strand::Strand;
    use rustymotif_utils::modtype::ModType;
    use std::fs::File;
    use std::io::Write;
    use rustymotif_utils::motif::Motif;
    use std::collections::HashMap;

    #[test]
    fn test_load_contig_bin() {
        let contig_bin_path = Path::new("test_data/contig_bin.tsv");
        let bin_contigs = load_contig_bin(contig_bin_path).expect("Failed to load contig bin");
        assert!(!bin_contigs.is_empty());
        assert!(bin_contigs.contains_key("bin1"));
        assert!(bin_contigs["bin1"].contains(&"contig1".to_string()));
        assert!(bin_contigs["bin1"].contains(&"contig2".to_string()));
    }
    #[test]
    fn test_load_bin_motifs() {
        let bin_motifs_path = Path::new("test_data/bin-motifs.tsv");
        let bin_motifs = load_bin_motifs(bin_motifs_path).expect("Failed to load bin motifs");
        assert!(!bin_motifs.is_empty());
        assert!(bin_motifs.contains_key("bin1"));
        let motifs = &bin_motifs["bin1"];
        assert_eq!(motifs.len(), 7);
    }   

    #[test]
    fn test_annotated_record_writer() {
        let file_path = Path::new("test_data/annotated_records.tsv");
        let mut writer = AnnotatedRecordWriter::new(file_path).expect("Failed to create writer");
        writer.write_header().expect("Failed to write header"); 
        let record = AnnotatedPileupRecord::new(
            PileupRecord {
                reference: "contig1".to_string(),
                position: 100,
                strand: Strand::Positive,
                mod_type: ModType::SixMA,
                n_mod: 2,
                n_valid_cov: 10,
                n_canonical: 8,
                n_diff: 2,
            },
            true,
            vec![Motif::new("GATC", "6mA", 1).unwrap()],
        );
        writer.write_record(&record).expect("Failed to write record");
        writer.flush().expect("Failed to sync file");

        // Check file contain two lines: header and record
        let content = fs::read_to_string(file_path).expect("Failed to read file");
        let lines: Vec<&str> = content.lines().collect();
        assert_eq!(lines.len(), 2, "File should contain exactly two lines: header and record");
        assert!(lines[1].contains("contig1"));
        assert!(lines[1].contains("6mA"));
        assert!(lines[1].contains("100"));
        assert!(lines[1].contains("true"));
        assert!(lines[1].contains("G6mATC")); // Check that motifs field is empty
    }
    #[test]
    fn test_annotated_record_writer_empty_motifs() {
        let file_path = Path::new("test_data/annotated_records_empty_motifs.tsv");
        let mut writer = AnnotatedRecordWriter::new(file_path).expect("Failed to create writer");
        writer.write_header().expect("Failed to write header");
        let record = AnnotatedPileupRecord::new(
            PileupRecord {
                reference: "contig1".to_string(),
                position: 100,  
                strand: Strand::Positive,
                mod_type: ModType::SixMA,
                n_mod: 2,
                n_valid_cov: 10,
                n_canonical: 8,
                n_diff: 2,
            },
            false,
            Vec::new(), 
        );
        writer.write_record(&record).expect("Failed to write record");
        writer.flush().expect("Failed to sync file");

        // Check file contain two lines: header and record
        let content = fs::read_to_string(file_path).expect("Failed to read file");
        let lines: Vec<&str> = content.lines().collect();
        assert_eq!(lines.len(), 2, "File should contain exactly two lines: header and record");
        assert!(lines[1].contains("contig1"));
        assert!(lines[1].contains("6mA"));
        assert!(lines[1].contains("100"));
        assert!(lines[1].contains("false"));
        assert!(lines[1].contains("")); // Check that motifs field is empty
    }
    #[test]
    fn test_annotated_record_writer_multiple_records() {
        let file_path = Path::new("test_data/annotated_records_multiple.tsv");
        let mut writer = AnnotatedRecordWriter::new(file_path).expect("Failed to create writer");
        writer.write_header().expect("Failed to write header");
        let records = vec![
            AnnotatedPileupRecord::new(
                PileupRecord {
                    reference: "contig1".to_string(),
                    position: 100,
                    strand: Strand::Positive,
                    mod_type: ModType::SixMA,
                    n_mod: 2,
                    n_valid_cov: 10,
                    n_canonical: 8,
                    n_diff: 2,
                },
                true,
                vec![Motif::new("GATC", "6mA", 1).unwrap()],
            ),
            AnnotatedPileupRecord::new(
                PileupRecord {
                    reference: "contig1".to_string(),
                    position: 200,
                    strand: Strand::Negative,
                    mod_type: ModType::FiveMC,
                    n_mod: 3,
                    n_valid_cov: 15,
                    n_canonical: 12,
                    n_diff: 3,  
                },
                false,
                vec![],
            ),
        ];
        writer.write_records_iter(records.iter()).expect("Failed to write records");
        writer.flush().expect("Failed to sync file");
        // Check file contain three lines: header and two records
        let content = fs::read_to_string(file_path).expect("Failed to read file");
        let lines: Vec<&str> = content.lines().collect();
        assert_eq!(lines.len(), 3, "File should contain exactly three lines: header and two records");
        assert!(lines[1].contains("contig1"));
        assert!(lines[1].contains("6mA"));
        assert!(lines[1].contains("100"));
        assert!(lines[1].contains("true"));
        assert!(lines[1].contains("G6mATC")); // Check that motifs
        assert!(lines[2].contains("contig1"));
        assert!(lines[2].contains("5mC"));
        assert!(lines[2].contains("200"));
        assert!(lines[2].contains("false"));
        assert!(lines[2].contains("")); // Check that motifs field is empty

    }
}
