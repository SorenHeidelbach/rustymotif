use crate::pileup::PileupChunk;
use crate::pileup::PileupRecord;
use ahash::{HashMap, HashMapExt, HashSet, HashSetExt};
use rustymotif_utils::{modtype::ModType, motif::Motif, pileup, strand::Strand};
use log::{debug, info};
use anyhow::{anyhow, Ok};


pub struct PileupChunkHashMap {
    pub contig: String,
    pub records: HashMap<(usize, Strand, ModType), PileupRecord>,
}

impl PileupChunkHashMap {
    pub fn new(contig: String) -> Self {
        Self {
            contig,
            records: HashMap::new(),
        }
    }

    pub fn from_chunk(chunk: PileupChunk) -> Self {
        let mut records = HashMap::new();
        for record in chunk.records {
            let key = (record.position, record.strand, record.mod_type);
            records.insert(key, record);
        }
        Self {
            contig: chunk.reference,
            records,
        }
    }

    pub fn add_record(&mut self, record: PileupRecord) {
        let key = (record.position, record.strand, record.mod_type);
        self.records.insert(key, record);
    }


    pub fn add_records(&mut self, mut records: PileupChunk) {
        assert_eq!(records.reference, "Reference mismatch");
        self.records.reserve(records.records.len());
        records.records.drain(..).for_each(|record| self.add_record(record));
    }



    pub fn get_record(&self, position: usize, strand: Strand, mod_type: ModType) -> Option<&PileupRecord> {
        self.records.get(&(position, strand, mod_type))
    }

    pub fn len(&self) -> usize {
        self.records.len()
    }

    pub fn retain_max_mod_per_pos_strand(&mut self) {
        // Step 1: Find the max n_mod per (pos, strand), ignoring ModType
        let mut max_mod_map: HashMap<(usize, Strand), u32> = HashMap::new();

        for ((pos, strand, _mod_type), record) in &self.records {
            let key = (*pos, *strand);
            max_mod_map
                .entry(key)
                .and_modify(|max| *max = (*max).max(record.n_mod))
                .or_insert(record.n_mod);
        }

        // Step 2: Retain only records with max n_mod at that (pos, strand)
        self.records
            .retain(|(pos, strand, _), record| {
                if let Some(&max_mod) = max_mod_map.get(&(*pos, *strand)) {
                    record.n_mod == max_mod
                } else {
                    false
                }
            });
    }

    pub fn remove_windows_with_n_high_methylation(
        &mut self,
        window_span: usize,
        fraction_threshold: f32,
        n_high_methylation: usize,
    ) {
        // Step 1: Group records by (ModType, Strand)
        let mut grouped: HashMap<(ModType, Strand), Vec<(usize, (usize, Strand, ModType), &PileupRecord)>> =
            HashMap::new();
        // debug!("Grouping records by (ModType, Strand)");
        for (key @ (pos, strand, mod_type), record) in &self.records {
            grouped
                .entry((*mod_type, *strand))
                .or_default()
                .push((*pos, *key, record));
        }

        let mut keys_to_remove = HashSet::new();

        // Step 2: For each group, do position-based sliding window
        // debug!("Processing grouped records for high methylation windows");
        for ((_mod_type, _strand), mut records) in grouped {
            // Sort by position
            records.sort_by_key(|(pos, _, _)| *pos);

            let n = records.len();
            for i in 0..n {
                let (start_pos, _, _) = records[i];

                // Extend window to include all records within window_span
                let mut j = i;
                while j < n && records[j].0 - start_pos <= window_span {
                    j += 1;
                }

                let window = &records[i..j];

                let high_frac_count = window.iter().filter(|(_, _, rec)| {
                    rec.n_valid_cov > 0 &&
                    (rec.n_mod as f32 / rec.n_valid_cov as f32) > fraction_threshold
                }).count();

                if high_frac_count >= n_high_methylation {
                    for (_, key, _) in window {
                        keys_to_remove.insert(*key);
                    }
                }
            }
        }

        // Step 3: Remove matching entries from self.records
        // debug!("Removing {} records with high methylation", keys_to_remove.len());
        self.records.retain(|key, _| !keys_to_remove.contains(key));
    }
}





#[cfg(test)]
mod tests {
    use super::*;
    use crate::pileup::PileupRecord;
    use rustymotif_utils::strand::Strand;
    use rustymotif_utils::modtype::ModType;

}



