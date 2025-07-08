use std::collections::HashSet;

use ahash::{HashMap, HashMapExt};
use anyhow::anyhow;
use anyhow::bail;
use anyhow::Result;
use bio::io::fasta;
use regex::Regex;
use utils::iupac::{self, IupacBase};
use utils::modtype::ModType;
use utils::strand::Strand;
use utils::motif::{Motif, ComplementMotif, MotifLike};

use utils::pileup::{PileupChunk, PileupRecord};
use serde::{Deserialize, Serialize};


#[cfg(test)]
mod tests {
    use super::*;
    use crate::motif::Motif;
    use crate::pileup::PileupRecord;
    use bio::bio_types::annot::contig;
    use utils::modtype::ModType;
    use utils::strand::Strand;

    #[test]
    fn test_contig_add_record() {
        let mut contig = Contig::new("test", "ACGT");
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
        let mut contig = Contig::new("test", "ACGT");
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
    fn test_contig_find_motif_indices() -> Result<()> {
        let contig = Contig::new("test", "ACGTACGTACGTACGT");
        let motif = Motif::new("ACGT", "6mA", 0).unwrap();
        let indices = contig.find_motif_indices(&motif)?.unwrap();
        assert_eq!(indices, vec![0, 4, 8, 12]);

        let contig = Contig::new(
            "test",
            "ACCCCGGAGGTCGTACGCCGGATCCGGTACCGGACGTACCGGTCGCCGGAT",
        );
        let motif = Motif::new("CCGGA", "6mA", 4).unwrap();
        let indices = contig.find_motif_indices(&motif)?;
        assert_eq!(indices, Some(vec![7, 21, 33, 49]));
        Ok(())
    }




    #[test]
    fn test_contig_extract_subsequences() {
        let contig = Contig::new("test", "ACGTACGTACGTACGTA");
        let indices = vec![2, 6, 10, 14];
        let subsequences = contig.extract_subsequences(&indices, 2).unwrap();
        assert_eq!(subsequences, vec!["ACGTA", "ACGTA", "ACGTA", "ACGTA"]);

        let contig = Contig::new("test", "ACGTACGTACGTACGT");
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
