use std::collections::HashSet;

use crate::motif::{Motif, MotifLike};
use crate::pileup::{PileupChunk, PileupRecord};
use crate::search::MotifNode;
use ahash::{HashMap, HashMapExt};
use anyhow::anyhow;
use anyhow::bail;
use anyhow::Result;
use bio::io::fasta;
use regex::Regex;
use rustymotif_utils::iupac::{self, IupacBase};
use rustymotif_utils::modtype::ModType;
use rustymotif_utils::strand::Strand;
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub enum MethylationLevel {
    High,
    Middle,
    Low,
}

#[derive(Debug, Clone, PartialEq)]
pub struct MethylationThresholds {
    pub middle: f32,
    pub high: f32,
}

impl MethylationThresholds {
    pub fn new(middle: f32, high: f32) -> Self {
        Self { middle, high }
    }

    pub fn default() -> Self {
        Self {
            middle: 0.2,
            high: 0.7,
        }
    }
}

#[derive(Debug, Clone)]
pub struct Contig {
    pub reference: String,
    pub sequence: String,
    pub records: HashMap<(usize, Strand, ModType), PileupRecord>,
    pub methylation_levels: HashMap<(MethylationLevel, Strand, ModType), Vec<usize>>,
    pub thresholds: MethylationThresholds,
}

impl Contig {
    pub fn new(reference: &str, sequence: &str, methylation_thresholds: MethylationThresholds) -> Self {
        Self {
            reference: reference.to_string(),
            sequence: sequence.to_string(),
            records: HashMap::new(),
            methylation_levels: HashMap::new(),
            thresholds: methylation_thresholds,
        }
    }


    pub fn add_record(&mut self, record: PileupRecord) {
        let key = (record.position, record.strand, record.mod_type);
        self.records.insert(key, record);
    }


    pub fn add_records(&mut self, mut records: PileupChunk) {
        assert_eq!(records.reference, self.reference, "Reference mismatch");
        self.records.reserve(records.records.len());
        records.records.drain(..).for_each(|record| self.add_record(record));
    }

    pub fn find_motif_indices(&self, motif: &Motif) -> Result<Option<Vec<usize>>> {
        let re = Regex::new(&motif.regex()?)?;
        let mut indices = Vec::new();
        re.find_iter(&self.sequence)
            .for_each(|m| indices.push(m.start() + motif.position as usize));
        Ok((!indices.is_empty()).then_some(indices))
    }

    pub fn find_complement_motif_indices(&self, motif: &Motif) -> Result<Option<Vec<usize>>> {
        let re = Regex::new(&motif.reverse_complement()?.regex()?)?;
        let mut indices = Vec::new();
        re.find_iter(&self.sequence)
            .for_each(|m| indices.push(m.start() + motif.position as usize));
        Ok((!indices.is_empty()).then_some(indices))
    }


    pub fn extract_subsequences(&self, indices: &[usize], flank_size: usize) -> Result<Vec<String>> {
        let max_index = self.sequence.len().saturating_sub(1);
        Ok(indices
            .iter()
            .filter_map(|&index| {
                if index >= flank_size && index + flank_size <= max_index {
                    Some(self.sequence[index - flank_size..=index + flank_size].to_owned())
                } else {
                    None
                }
            })
            .collect())
    }

    pub fn get_by_mod_type(&self, mod_type: ModType) -> Vec<&PileupRecord> {
        self.records
            .iter()
            .filter_map(|(key, record)| {
                if key.2 == mod_type {
                    Some(record)
                } else {
                    None
                }
            })
            .collect()
    }

    pub fn remove_record_at_motif_indices(
        &mut self,
        motif: &Motif,
    ) -> Result<Vec<PileupRecord>> {
        let indices = self.find_motif_indices(motif)?;
        if let Some(indices) = indices {
            let mut removed_records = Vec::new();
            for index in indices {
                if let Some(record) = self.records.remove(&(index, Strand::Positive, motif.mod_type)) {
                    removed_records.push(record);
                }
            }
            Ok(removed_records)
        } else {
            bail!("No motif found");
        }
    }
    
    pub fn get_records_for_motif(&self, motif: &Motif) -> Result<Vec<&PileupRecord>> {
        let mut records = vec![];
        let sequence_chars: Vec<char> = self.sequence.chars().collect();
        let motif_len = motif.sequence.len();
        
        // Search for motif matches considering IUPAC ambiguity codes
        for i in 0..=(sequence_chars.len().saturating_sub(motif_len)) {
            let mut matches = true;
            
            // Check if the sequence at this position matches the motif
            for (j, motif_base) in motif.sequence.iter().enumerate() {
                if let Ok(seq_base) = IupacBase::from_char(sequence_chars[i + j]) {
                    if !motif_base.contain_other_iupac_base(seq_base) {
                        matches = false;
                        break;
                    }
                } else {
                    matches = false;
                    break;
                }
            }
            
            if matches {
                let motif_pos = i + motif.position as usize;
                // Look for records at this position for the motif's modification type
                let key = (motif_pos, Strand::Positive, motif.mod_type);
                if let Some(record) = self.records.get(&key) {
                    records.push(record);
                }
                // Also check negative strand
                let key = (motif_pos, Strand::Negative, motif.mod_type);
                if let Some(record) = self.records.get(&key) {
                    records.push(record);
                }
            }
        }
        Ok(records)
    }
    pub fn remove_records_for_motif(&mut self, motif: &Motif) -> Result<()> {
        let records = self.get_records_for_motif(motif)?;
        let keys = records
            .iter()
            .map(|record| (record.position, record.strand, record.mod_type))
            .collect::<HashSet<_>>();
        drop(records); // Release the borrow on records before modifying self.records
        for key in keys {
            self.records.remove(&key);
        }
        Ok(())
    }
    
    pub fn classify_methylation_level(&self, n_mod: f32, n_valid_cov: f32) -> MethylationLevel {
        if n_mod == 0.0  {
            return MethylationLevel::Low;
        }
        if n_valid_cov == 0.0 {
            return MethylationLevel::Low;
        }
        match n_mod / n_valid_cov {
            n if n >= self.thresholds.high => MethylationLevel::High,
            n if n >= self.thresholds.middle => MethylationLevel::Middle,
            _ => MethylationLevel::Low,
        }
    }
    
    pub fn populate_methylation_levels(&mut self) {
        self.methylation_levels.clear();
        for record in self.records.values() {
            let methylation_level = self.classify_methylation_level(record.n_mod as f32, record.n_valid_cov as f32);
            self.methylation_levels
                .entry((methylation_level, record.strand, record.mod_type))
                .or_default()
                .push(record.position);
        }
    }
    
    pub fn methylated_motif_positions(
        &self,
        motif: &Motif,
        methylation_level: MethylationLevel,
        strand: Strand,
    ) -> Result<Option<Vec<usize>>> {
        let Some(motif_indices) = self.find_motif_indices(motif)? else {
            return Ok(None);
        };
    
        let Some(methylated_indices) = self.methylation_levels.get(&(methylation_level, strand, motif.mod_type)) else {
            return Ok(None);
        };
        let motif_set: HashSet<_> = motif_indices.iter().copied().collect();
        let methyl_set: HashSet<_> = methylated_indices.iter().copied().collect();

        let intersection = motif_set
            .intersection(&methyl_set)
            .copied()
            .collect::<Vec<_>>();

        Ok(Some(intersection))
    }
}
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct EqualLengthDNASet {
    pub sequences: Vec<String>,
    pub length: usize,
}

impl EqualLengthDNASet {
    pub fn new_empty(length: usize) -> Result<Self> {
        Ok(Self {
            sequences: Vec::new(),
            length: length,
        })
    }

    pub fn new(sequences: Vec<String>) -> Result<Self> {
        if sequences.is_empty() {
            bail!("No sequences provided");
        }
        let length = sequences[0].len();
        if sequences.iter().all(|s| s.len() == length) {
            Ok(Self { sequences, length })
        } else {
            bail!("Sequences are not all the same length");
        }
    }

    pub fn len(&self) -> usize {
        self.sequences.len()
    }

    pub fn add_other(&mut self, other: &Self) -> Result<()> {
        if self.length != other.length {
            bail!("Mismatch in sequence length of the two sets");
        }
        self.sequences.extend_from_slice(&other.sequences);
        Ok(())
    }

    pub fn reverse_complement(&self) -> Self {
        let sequences = self
            .sequences
            .iter()
            .map(|seq| {
                seq.chars()
                    .rev()
                    .map(|c| IupacBase::from_char(c).unwrap().complement().to_string())
                    .collect()
            })
            .collect();
        Self {
            sequences,
            length: self.length,
        }
    }

    pub fn pssm(&self, pseudocount: f64) -> PSSM {
        let n = self.len() as f64;
        let length = self.length;

        // We'll accumulate counts in a 4 x length matrix
        let mut counts = vec![vec![0.0; length]; 4];

        // Tally
        for seq in &self.sequences {
            for (i, c) in seq.as_str().chars().enumerate() {
                match c {
                    'A' => counts[0][i] += 1.0,
                    'C' => counts[1][i] += 1.0,
                    'G' => counts[2][i] += 1.0,
                    'T' => counts[3][i] += 1.0,
                    _ => {
                        // If you handle 'N' or IUPAC codes, do that here
                    }
                }
            }
        }

        // Convert counts to frequencies, adding pseudocount
        let mut pssm = vec![vec![0.0; length]; 4];
        for b in 0..4 {
            for i in 0..length {
                pssm[b][i] = (counts[b][i] + pseudocount) / (n + 4.0 * pseudocount);
            }
        }

        PSSM::new(pssm)
    }

    // Extract sequence with alignment to the motif
    // The motif methylation position should be contered in the seqnece for it to match
    pub fn get_motif_matching_sequences(&self, motif: &Motif) -> Result<EqualLengthDNASet> {
        let seq_len = self.length;
        let anchor_idx = self.length / 2;
        let motif_len = motif.sequence.len();
        let motif_start = anchor_idx.saturating_sub(motif.position as usize);
        let motif_end = motif_start + motif_len;
        if motif_end > seq_len {
            anyhow::bail!("Motif is too long for the sequence");
        }

        let filtered_seq: Vec<String> = self
            .sequences
            .iter()
            .filter(|seq| {
                let seq_slice = &seq[motif_start..motif_end];
                seq_slice
                    .chars()
                    .zip(&motif.sequence)
                    .all(|(seq_base, motif_base)| {
                        motif_base.contain_other_iupac_base(
                            iupac::IupacBase::from_char(seq_base).unwrap(),
                        )
                    })
            })
            .cloned()
            .collect();
        if filtered_seq.is_empty() {
            bail!("No matching sequences found");
        }
        Ok(Self::new(filtered_seq)?)
    }

    pub fn remove_motif_matching_sequences(&self, motif: &Motif) -> Result<EqualLengthDNASet> {
        let seq_len = self.length;
        let anchor_idx = self.length / 2;
        let motif_len = motif.sequence.len();
        let motif_start = anchor_idx.saturating_sub(motif.position as usize);
        let motif_end = motif_start + motif_len;
        if motif_end > seq_len {
            anyhow::bail!("Motif is too long for the sequence");
        }

        let filtered_seq: Vec<String> = self
            .sequences
            .iter()
            .filter(|seq| {
                let seq_slice = &seq[motif_start..motif_end];
                !seq_slice
                    .chars()
                    .zip(&motif.sequence)
                    .all(|(seq_base, motif_base)| {
                        motif_base.contain_other_iupac_base(
                            iupac::IupacBase::from_char(seq_base).unwrap(),
                        )
                    })
            })
            .cloned()
            .collect();
        if filtered_seq.is_empty() {
            bail!(
                "No sequences remaining after filtering with {:?} motif",
                motif.sequence
            );
        }
        Ok(Self::new(filtered_seq)?)
    }
}

// Position Specific Scoring Matrix
// A PSSM is a 4 x length matrix of probabilities
// where each row corresponds to a base (A, C, G, T)
// and each column corresponds to a position in a sequence
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct PSSM {
    pssm: Vec<Vec<f64>>,
    length: usize,
}

impl PSSM {
    pub fn new(pssm: Vec<Vec<f64>>) -> Self {
        let length = pssm[0].len();
        Self { pssm, length }
    }

    pub fn length(&self) -> usize {
        self.length
    }

    pub fn kl_divergence(&self, other: &PSSM) -> Result<Vec<f64>> {
        if self.length != other.length {
            bail!("Length mismatch");
        }
        let mut kl_div = Vec::new();
        let epsilon = 1e-6;
        for i in 0..self.length {
            let mut kl_div_i = 0.0;
            for b in 0..4 {
                let p = self.pssm[b][i] + epsilon;
                let q = other.pssm[b][i] + epsilon;
                kl_div_i += p * (p / q).ln();
            }
            kl_div.push(kl_div_i);
        }
        Ok(kl_div)
    }

    pub fn bases_above_freq(&self, position: u8, min_freq: f64) -> Vec<IupacBase> {
        let possible_bases = vec![IupacBase::A, IupacBase::C, IupacBase::G, IupacBase::T];
        let mut bases = Vec::new();
        for i in 0..4 {
            let base = possible_bases[i];
            if self.pssm[i][position as usize] > min_freq {
                bases.push(base);
            }
        }
        bases
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::motif::Motif;
    use crate::pileup::PileupRecord;
    use bio::bio_types::annot::contig;
    use rustymotif_utils::modtype::ModType;
    use rustymotif_utils::strand::Strand;

    #[test]
    fn test_contig_add_record() {
        let mut contig = Contig::new("test", "ACGT", MethylationThresholds::default());
        let record = PileupRecord {
            reference: "test".to_string(),
            position: 0,
            strand: Strand::Positive,
            mod_type: ModType::SixMA,
            n_mod: 1,
            n_valid_cov: 1,
            n_canonical: 1,
            n_diff: 0,
        };
        contig.add_record(record.clone());
        assert_eq!(contig.records.len(), 1);
        assert_eq!(
            contig.records.get(&(0, Strand::Positive, ModType::SixMA)),
            Some(&record)
        );
    }

    #[test]
    fn test_contig_add_records() {
        let mut contig = Contig::new("test", "ACGT", MethylationThresholds::default());
        let records = PileupChunk {
            reference: "test".to_string(),
            records: vec![
                PileupRecord {
                    reference: "test".to_string(),
                    position: 0,
                    strand: Strand::Positive,
                    mod_type: ModType::SixMA,
                    n_mod: 1,
                    n_valid_cov: 1,
                    n_canonical: 1,
                    n_diff: 0,
                },
                PileupRecord {
                    reference: "test".to_string(),
                    position: 1,
                    strand: Strand::Positive,
                    mod_type: ModType::SixMA,
                    n_mod: 1,
                    n_valid_cov: 1,
                    n_canonical: 1,
                    n_diff: 0,
                },
            ],
        };
        contig.add_records(records.clone());
        assert_eq!(contig.records.len(), 2);
        assert_eq!(
            contig.records.get(&(0, Strand::Positive, ModType::SixMA)),
            Some(&records.records[0])
        );
        assert_eq!(
            contig.records.get(&(1, Strand::Positive, ModType::SixMA)),
            Some(&records.records[1])
        );
    }

    #[test]
    fn test_contig_populate_methylation_levels() {
        let mut contig = Contig::new("test", "ACGT", MethylationThresholds::default());
        let records = PileupChunk {
            reference: "test".to_string(),
            records: vec![
                PileupRecord {
                    reference: "test".to_string(),
                    position: 0,
                    strand: Strand::Positive,
                    mod_type: ModType::SixMA,
                    n_mod: 1,
                    n_valid_cov: 1,
                    n_canonical: 1,
                    n_diff: 0,
                },
                PileupRecord {
                    reference: "test".to_string(),
                    position: 1,
                    strand: Strand::Positive,
                    mod_type: ModType::SixMA,
                    n_mod: 0,
                    n_valid_cov: 1,
                    n_canonical: 1,
                    n_diff: 0,
                },
            ],
        };
        contig.add_records(records.clone());
        contig.populate_methylation_levels();
        assert_eq!(
            contig.methylation_levels.get(&(MethylationLevel::High, Strand::Positive, ModType::SixMA)),
            Some(&vec![0])
        );
        assert_eq!(
            contig.methylation_levels.get(&(MethylationLevel::Low, Strand::Positive, ModType::SixMA)),
            Some(&vec![1])
        );
    }

    #[test]
    fn test_contig_find_motif_indices() -> Result<()> {
        let contig = Contig::new("test", "ACGTACGTACGTACGT", MethylationThresholds::default());
        let motif = Motif::new("ACGT", "6mA", 0).unwrap();
        let indices = contig.find_motif_indices(&motif)?.unwrap();
        assert_eq!(indices, vec![0, 4, 8, 12]);

        let contig = Contig::new(
            "test",
            "ACCCCGGAGGTCGTACGCCGGATCCGGTACCGGACGTACCGGTCGCCGGAT",
            MethylationThresholds::default(),
        );
        let motif = Motif::new("CCGGA", "6mA", 4).unwrap();
        let indices = contig.find_motif_indices(&motif)?;
        assert_eq!(indices, Some(vec![7, 21, 33, 49]));
        Ok(())
    }




    #[test]
    fn test_contig_extract_subsequences() {
        let contig = Contig::new("test", "ACGTACGTACGTACGTA", MethylationThresholds::default());
        let indices = vec![2, 6, 10, 14];
        let subsequences = contig.extract_subsequences(&indices, 2).unwrap();
        assert_eq!(subsequences, vec!["ACGTA", "ACGTA", "ACGTA", "ACGTA"]);

        let contig = Contig::new("test", "ACGTACGTACGTACGT", MethylationThresholds::default());
        let indices = vec![0, 4, 8, 12];
        let subsequences = contig.extract_subsequences(&indices, 3).unwrap();
        assert_eq!(subsequences, vec!["CGTACGT", "CGTACGT", "CGTACGT"]);
    }

    #[test]
    fn test_equal_length_dna_set_new() {
        let sequences = vec!["ACGT".to_string(), "ACGT".to_string()];
        let dna_set = EqualLengthDNASet::new(sequences);
        assert!(dna_set.is_ok());

        let sequences = vec!["ACGT".to_string(), "ACG".to_string()];
        let dna_set = EqualLengthDNASet::new(sequences);
        assert!(dna_set.is_err());
    }

    #[test]
    fn test_equal_length_dna_set_reverse_complement() {
        let sequences = vec!["ACGT".to_string(), "ACGT".to_string()];
        let dna_set = EqualLengthDNASet::new(sequences).unwrap();
        let rc_dna_set = dna_set.reverse_complement();
        assert_eq!(
            rc_dna_set.sequences,
            vec!["ACGT".to_string(), "ACGT".to_string()]
        );

        let sequences = vec!["TACG".to_string(), "ACGT".to_string()];
        let dna_set = EqualLengthDNASet::new(sequences).unwrap();
        let rc_dna_set = dna_set.reverse_complement();
        assert_eq!(
            rc_dna_set.sequences,
            vec!["CGTA".to_string(), "ACGT".to_string()]
        );
    }

    #[test]
    fn test_equal_length_dna_set_empty() {
        let dna_set = EqualLengthDNASet::new_empty(4);
        assert!(dna_set.is_ok());
    }

    #[test]
    fn test_equal_length_dna_set_get_motif_matching_sequences() {
        let sequences = vec![
            "ACGTACGTACGTACGG".to_string(),
            "AAAAWAAAAAAAAAAA".to_string(),
            "ACGTACGTACGTACGT".to_string(),
            "ACGTACGTACGTACGT".to_string(),
        ];
        let dna_set = EqualLengthDNASet::new(sequences).unwrap();
        let motif = Motif::new("ACGT", "6mA", 0).unwrap();
        let matching_sequences = dna_set.get_motif_matching_sequences(&motif).unwrap();
        assert_eq!(
            matching_sequences,
            EqualLengthDNASet::new(vec![
                "ACGTACGTACGTACGG".to_string(),
                "ACGTACGTACGTACGT".to_string(),
                "ACGTACGTACGTACGT".to_string(),
            ])
            .unwrap()
        );

        let sequences = vec![
            "AAAAWAAAAAAAAAAA".to_string(),
            "AAAAWAAAAAAAAAAA".to_string(),
            "AAAAWAAAAAAAAAAA".to_string(),
            "AAAAWAAAAAAAAAAA".to_string(),
        ];
        let dna_set = EqualLengthDNASet::new(sequences).unwrap();
        let motif = Motif::new("ACGT", "6mA", 0).unwrap();
        let matching_sequences = dna_set.get_motif_matching_sequences(&motif);
        assert!(matching_sequences.is_err());
    }

    #[test]
    fn test_equal_length_dna_set_pssm() {
        let sequences = vec![
            "ACGT".to_string(),
            "ACGT".to_string(),
            "ACGT".to_string(),
            "ACGT".to_string(),
        ];
        let dna_set = EqualLengthDNASet::new(sequences).unwrap();
        let pssm = dna_set.pssm(0.0);
        assert_eq!(
            pssm,
            PSSM::new(vec![
                vec![1.0, 0.0, 0.0, 0.0],
                vec![0.0, 1.0, 0.0, 0.0],
                vec![0.0, 0.0, 1.0, 0.0],
                vec![0.0, 0.0, 0.0, 1.0],
            ])
        );

        let sequences = vec![
            "TACG".to_string(),
            "ACGT".to_string(),
            "ACGT".to_string(),
            "ACGT".to_string(),
        ];
        let dna_set = EqualLengthDNASet::new(sequences).unwrap();
        let pssm = dna_set.pssm(0.0);
        assert_eq!(
            pssm,
            PSSM::new(vec![
                vec![0.75, 0.25, 0.0, 0.0],
                vec![0.0, 0.75, 0.25, 0.0],
                vec![0.0, 0.0, 0.75, 0.25],
                vec![0.25, 0.0, 0.0, 0.75],
            ])
        );
        let sequences = vec![
            "TACGG".to_string(),
            "ACGTG".to_string(),
            "ACGTG".to_string(),
            "ACGTG".to_string(),
        ];
        let dna_set = EqualLengthDNASet::new(sequences).unwrap();
        let pssm = dna_set.pssm(0.0);
        assert_eq!(
            pssm,
            PSSM::new(vec![
                vec![0.75, 0.25, 0.0, 0.0, 0.0],
                vec![0.0, 0.75, 0.25, 0.0, 0.0],
                vec![0.0, 0.0, 0.75, 0.25, 1.0],
                vec![0.25, 0.0, 0.0, 0.75, 0.0],
            ])
        );
    }

    #[test]
    fn test_pssm_kl_divergence() {
        let pssm1 = PSSM::new(vec![
            vec![0.75, 0.25, 0.0, 0.0],
            vec![0.0, 0.75, 0.25, 0.0],
            vec![0.0, 0.0, 0.75, 0.25],
            vec![0.25, 0.0, 0.0, 0.75],
        ]);
        let pssm2 = PSSM::new(vec![
            vec![0.75, 0.25, 0.0, 0.0],
            vec![0.0, 0.75, 0.25, 0.0],
            vec![0.0, 0.0, 0.75, 0.25],
            vec![0.25, 0.0, 0.0, 0.75],
        ]);
        assert_eq!(
            pssm1.kl_divergence(&pssm2).unwrap(),
            vec![0.0, 0.0, 0.0, 0.0]
        );

        let pssm1 = PSSM::new(vec![
            vec![0.75, 0.25, 0.0, 0.0],
            vec![0.0, 0.75, 0.25, 0.0],
            vec![0.0, 0.0, 0.75, 0.25],
            vec![0.25, 0.0, 0.0, 0.75],
        ]);
        let pssm2 = PSSM::new(vec![
            vec![1.0, 0.0, 0.0, 0.0],
            vec![0.0, 1.0, 0.0, 0.0],
            vec![0.0, 0.0, 1.0, 0.0],
            vec![0.0, 0.0, 0.0, 1.0],
        ]);
        assert_eq!(
            pssm1.kl_divergence(&pssm2).unwrap(),
            vec![
                2.891555886408426,
                2.891555886408426,
                2.891555886408426,
                2.891555886408426
            ]
        );
    }
}
