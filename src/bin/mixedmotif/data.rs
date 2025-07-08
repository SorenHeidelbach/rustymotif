use crate::pileup::PileupChunk;
use crate::sequence::Contig;
use crate::sequence::MethylationThresholds;
use ahash::{HashMap, HashMapExt};
use log::debug;
pub struct GenomeWorkSpaceBuilder {
    pub contigs: HashMap<String, Contig>,
}

impl GenomeWorkSpaceBuilder {
    pub fn new() -> Self {
        Self {
            contigs: HashMap::new(),
        }
    }

    pub fn add_contig(&mut self, contig: Contig) {
        self.contigs.insert(contig.reference.clone(), contig);
    }

    pub fn push_records(&mut self, records: PileupChunk) {
        let contig = self
            .contigs
            .get_mut(records.reference.as_str())
            .expect("Could not find contig");
        contig.add_records(records);
        debug!("Added {} records to contig", contig.records.len());
    }

    pub fn build(self) -> GenomeWorkspace {
        GenomeWorkspace {
            contigs: self.contigs,
        }
    }
}

pub struct GenomeWorkspace {
    pub contigs: HashMap<String, Contig>,
}


#[cfg(test)]
mod tests {
    use super::*;
    use crate::pileup::PileupRecord;
    use crate::sequence::{Contig};
    use rustymotif_utils::strand::Strand;
    use rustymotif_utils::modtype::ModType;

    #[test]
    fn test_genome_workspace_builder() {
        let mut builder = GenomeWorkSpaceBuilder::new();
        let contig = Contig::new("contig_1", "AACG", MethylationThresholds::new(0.1, 0.5));
        builder.add_contig(contig.clone());
        let records = PileupChunk {
            reference: "contig_1".to_string(),
            records: vec![
                PileupRecord {
                    reference: "contig_1".to_string(),
                    position: 0,
                    strand: Strand::Positive,
                    mod_type: ModType::SixMA,
                    n_mod: 1,
                    n_valid_cov: 1,
                    n_canonical: 1,
                    n_diff: 0,
                },
                PileupRecord {
                    reference: "contig_1".to_string(),
                    position: 1,
                    strand: Strand::Positive,
                    mod_type: ModType::SixMA,
                    n_mod: 0,
                    n_valid_cov: 1,
                    n_canonical: 1,
                    n_diff: 0,
                },
                PileupRecord {
                    reference: "contig_1".to_string(),
                    position: 2,
                    strand: Strand::Positive,
                    mod_type: ModType::FiveMC,
                    n_mod: 50,
                    n_valid_cov: 100,
                    n_canonical: 5,
                    n_diff: 0,
                },
                PileupRecord {
                    reference: "contig_1".to_string(),
                    position: 3,
                    strand: Strand::Positive,
                    mod_type: ModType::FiveMC,
                    n_mod: 50,
                    n_valid_cov: 100,
                    n_canonical: 5,
                    n_diff: 0,
                },
            ],
        };
        builder.push_records(records);
        let workspace = builder.build();
        assert_eq!(workspace.contigs.len(), 1);
        let contig = workspace.contigs.get("contig_1").unwrap();
        assert_eq!(contig.records.len(), 4);
    }
}