use crate::pileup::PileupChunk;
use crate::sequence::Contig;
use crate::pileup::PileupRecord;
use ahash::{HashMap, HashMapExt, HashSet};
use rustymotif_utils::{modtype::ModType, motif::Motif, pileup, strand::Strand};
use log::{debug, info};
use anyhow::{anyhow, Ok};

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

pub struct AnnotatedPileupRecord {
    pub pileup_record: PileupRecord,
    pub in_motif: bool,
    pub motifs: Vec<Motif>,
}

impl AnnotatedPileupRecord {
    pub fn new(
        pileup_record: PileupRecord,
        in_motif: bool,
        motifs: Vec<Motif>,
    ) -> Self {
        Self {
            pileup_record: pileup_record,
            in_motif,
            motifs,
        }
    }
}

pub struct AnnotatedPileupSummary {
    pub summary: HashMap<(Strand, ModType), Vec<String>>,
    pub contig_id: String,
}

impl AnnotatedPileupSummary {
    pub fn new() -> Self {
        Self {
            summary: HashMap::new(),
            contig_id: String::new(),
        }
    }

    pub fn summarise_records(
        &mut self,
        annotated_records: &[AnnotatedPileupRecord],
        min_cov: u32,
        threshold: f32,
    ) -> Result<(), anyhow::Error> {
        let mod_types: HashSet<ModType> = annotated_records
            .iter()
            .map(|r| r.pileup_record.mod_type)
            .collect();
        let contigs = annotated_records
            .iter()
            .map(|r| r.pileup_record.reference.clone())
            .collect::<HashSet<String>>();
        if contigs.len() != 1 {
            panic!("Annotated records must belong to a single contig");
        }
        self.contig_id = contigs.into_iter().next().unwrap();
        for mod_type in mod_types.into_iter() {
            for strand in [Strand::Positive, Strand::Negative] {
                let filtered_annotated_records = annotated_records.iter()
                    .filter(|r| r.pileup_record.mod_type == mod_type)
                    .filter(|r| r.pileup_record.strand == strand)
                    .filter(|r| r.pileup_record.n_valid_cov > min_cov)
                    .collect::<Vec<&AnnotatedPileupRecord>>();

                let filtered_annotated_records_in_motif = filtered_annotated_records.iter()
                    .filter(|r| r.in_motif)
                    .cloned()
                    .collect::<Vec<&AnnotatedPileupRecord>>();
                let filtered_annotated_records_not_in_motif = filtered_annotated_records.iter()
                    .filter(|r| !r.in_motif)
                    .cloned()
                    .collect::<Vec<&AnnotatedPileupRecord>>();

                // Extract unique motifs from the motifs in records
                let motifs_set: HashSet<Motif> = filtered_annotated_records_in_motif.iter()
                    .flat_map(|r| r.motifs.iter())
                    .cloned()
                    .collect();
                let n_motifs = motifs_set.len();

                let n_records = filtered_annotated_records.len();
                let n_records_in_motif = filtered_annotated_records_in_motif.iter()
                    .count();
                let n_records_not_in_motif = n_records - n_records_in_motif;
                let n_mod_in_motif = filtered_annotated_records_in_motif.iter()
                    .map(|r| r.pileup_record.n_mod)
                    .sum::<u32>();
                let n_canonical_in_motif = filtered_annotated_records_in_motif.iter()
                    .map(|r| r.pileup_record.n_canonical)
                    .sum::<u32>();
                let n_valid_cov_in_motif = filtered_annotated_records_in_motif.iter()
                    .map(|r| r.pileup_record.n_valid_cov)
                    .sum::<u32>();
                let n_diff_in_motif = filtered_annotated_records_in_motif.iter()
                    .map(|r| r.pileup_record.n_diff)
                    .sum::<u32>();
                let n_gt_threshold_in_motif = filtered_annotated_records_in_motif.iter()
                    .filter(|r| (r.pileup_record.n_mod as f32 / r.pileup_record.n_valid_cov as f32) > threshold)
                    .count();
                let n_mod_not_in_motif = filtered_annotated_records_not_in_motif.iter()
                    .filter(|r| !r.in_motif)
                    .map(|r| r.pileup_record.n_mod)
                    .sum::<u32>();
                let n_canonical_not_in_motif = filtered_annotated_records_not_in_motif.iter()
                    .map(|r| r.pileup_record.n_canonical)
                    .sum::<u32>();
                let n_valid_cov_not_in_motif = filtered_annotated_records_not_in_motif.iter()
                    .map(|r| r.pileup_record.n_valid_cov)
                    .sum::<u32>();
                let n_diff_not_in_motif = filtered_annotated_records_not_in_motif.iter()
                    .map(|r| r.pileup_record.n_diff)
                    .sum::<u32>();
                let n_gt_threshold_not_in_motif = filtered_annotated_records_not_in_motif.iter()
                    .filter(|r| (r.pileup_record.n_mod as f32 / r.pileup_record.n_valid_cov as f32) > threshold)
                    .count();
                let n_records_below_cov = annotated_records.iter()
                    .filter(|r| r.pileup_record.mod_type == mod_type)
                    .filter(|r| r.pileup_record.strand == strand)
                    .filter(|r| r.pileup_record.n_valid_cov <= min_cov)
                    .count();
                let mut summary = HashMap::new();
                summary.insert(
                    (strand, mod_type),
                    vec![
                        self.contig_id.clone(),
                        mod_type.to_string().to_string(),
                        strand.to_string(),
                        motifs_set.iter()
                            .map(|m| m.as_pretty_string())
                            .collect::<Vec<String>>()
                            .join(","),
                        n_motifs.to_string(),
                        n_records.to_string(),
                        n_records_below_cov.to_string(),
                        n_records_in_motif.to_string(),
                        n_records_not_in_motif.to_string(),
                        n_mod_in_motif.to_string(),
                        n_canonical_in_motif.to_string(),
                        n_valid_cov_in_motif.to_string(),
                        n_diff_in_motif.to_string(),
                        n_gt_threshold_in_motif.to_string(),
                        n_mod_not_in_motif.to_string(),
                        n_canonical_not_in_motif.to_string(),
                        n_valid_cov_not_in_motif.to_string(),
                        n_diff_not_in_motif.to_string(),
                        n_gt_threshold_not_in_motif.to_string(),
                    ],
                );
                self.summary.extend(summary);
            }
        }
        Ok(())
    }
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
        let contig = Contig::new("contig_1", "AACG");
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

