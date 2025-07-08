use crate::pileup::PileupChunk;
use crate::pileup::PileupRecord;
use crate::{modtype::ModType, motif::Motif, pileup, strand::Strand};
use ahash::{HashMap, HashMapExt, HashSet, HashSetExt};
use log::{debug, info, warn};
use anyhow::{anyhow};
use anyhow::Result;
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList};

#[pyclass]
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


    pub fn add_records(&mut self, mut records: PileupChunk) -> Result<()> {
        if records.reference != self.contig {
            return Err(anyhow!("Reference mismatch: expected {}, got {}", self.contig, records.reference));
        }
        self.records.reserve(records.records.len());
        records.records.drain(..).for_each(|record| self.add_record(record));
        Ok(())
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

#[pymodule]
fn rustymotif_utils(_py: Python, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PileupChunkHashMap>()?;
    m.add_function(wrap_pyfunction!(get_contig_pileup_data, m)?)?;
    Ok(())
}

#[pyfunction]
pub fn get_contig_pileup_data(
    py: Python<'_>,
    contig_id: &str,
    pileup_path: String,
    threads: usize,
    min_cov: u32,
    window_size: usize,
    n_high_methylation: usize,
    fraction_threshold: f32,
) -> PyResult<Py<PyDict>> {
    let mut chunk = load_contig(contig_id, pileup_path, threads)
        .map_err(|e| PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!("Failed to load contig: {}", e)))?;
    filter_chunk(&mut chunk, min_cov, window_size, n_high_methylation, fraction_threshold)
        .map_err(|e| PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!("Failed to filter chunk: {}", e)))?;

    // return polars friendly format
    // Create column-wise data    // Flatten records across all contigs
    let records: Vec<&PileupRecord> = chunk.records.values().collect();
    let dict = PyDict::new(py);
    dict.set_item("reference", records.iter().map(|r| r.reference.clone()).collect::<Vec<_>>())?;
    dict.set_item("position", records.iter().map(|r| r.position).collect::<Vec<_>>())?;
    dict.set_item("strand", records.iter().map(|r| r.strand.to_string()).collect::<Vec<_>>())?;
    dict.set_item("mod_type", records.iter().map(|r| r.mod_type.to_string()).collect::<Vec<_>>())?;
    dict.set_item("n_mod", records.iter().map(|r| r.n_mod).collect::<Vec<_>>())?;
    dict.set_item("n_valid_cov", records.iter().map(|r| r.n_valid_cov).collect::<Vec<_>>())?;
    dict.set_item("n_canonical", records.iter().map(|r| r.n_canonical).collect::<Vec<_>>())?;
    dict.set_item("n_diff", records.iter().map(|r| r.n_diff).collect::<Vec<_>>())?;

    Ok(dict.into())
}

pub fn load_contig(
    contig_id: &str,
    pileup_path: String,
    threads: usize,
) -> Result<PileupChunkHashMap> {
    let mut pileup_reader = pileup::PileupTabixReader::new(pileup_path, 5, threads).unwrap();
    info!("Loading contig {}", contig_id);
    pileup_reader.fetch_contig(contig_id).unwrap();
    match pileup_reader.read_all_from_contig() {
        Some(chunk_vec) => {
            let chunk = PileupChunkHashMap::from_chunk(chunk_vec);
            info!("Loaded {} records from contig {}", chunk.len(), contig_id);
            Ok(chunk)
        },
        None => {
            warn!("No records found for contig {}", contig_id);
            Err(anyhow!("No records found for contig {}", contig_id))
        }
    }

}

pub fn filter_chunk(
    chunk: &mut PileupChunkHashMap,
    min_cov: u32,
    window_size: usize,
    n_high_methylation: usize,
    fraction_threshold: f32,
) -> Result<()> {
    let original_len = chunk.len();
    info!("Remaining: {:5.1}%. Read {} records", 100.0*(chunk.records.len() as f64 /original_len as f64), chunk.records.len());
    chunk.retain_max_mod_per_pos_strand();
    if chunk.records.is_empty() {
        warn!("No records remaining after retaining max mod type per position and strand");
        return Err(anyhow!("No records remaining after retaining max mod type per position and strand"));
    }
    info!("Remaining: {:5.1}%. Retained only maximum mod type at all positions. {} records Remaining", 100.0*(chunk.records.len() as f64 /original_len as f64), chunk.len());
    chunk.records.retain(|_key, record| {
        record.n_valid_cov > min_cov
    });
    if chunk.records.is_empty() {
        warn!("No records remaining after retaining valid coverage > 5");
        return Err(anyhow!("No records remaining after retaining max mod type per position and strand"));
    }
    info!("Remaining: {:5.1}%. Retained only sufficient coverage positiont. Min cov: {}.: {} records Remaining", 100.0*(chunk.records.len() as f64 /original_len as f64), 5, chunk.len());
    chunk.remove_windows_with_n_high_methylation(
        window_size,
        fraction_threshold,
        n_high_methylation,
    );
    info!("Remaining: {:5.1}%. Removed windows with high methylation density. {} records Remaining",100.0*(chunk.records.len() as f64 /original_len as f64), chunk.len());
    info!("-- Finished processing");
    Ok(())
}




#[cfg(test)]
mod tests {
    use super::*;
    use crate::{modtype::ModType, strand::Strand};

    fn dummy_record(pos: usize, strand: Strand, mod_type: ModType, n_mod: u32, n_valid_cov: u32) -> PileupRecord {
        PileupRecord {
            reference: "contig1".to_string(),
            position: pos,
            strand,
            mod_type,
            n_mod,
            n_valid_cov,
            n_canonical: 0,
            n_diff: 0,
        }
    }

    #[test]
    fn test_new_empty_hashmap() {
        let map = PileupChunkHashMap::new("contig1".into());
        assert_eq!(map.len(), 0);
        assert_eq!(map.contig, "contig1");
    }

    #[test]
    fn test_add_record_and_get() {
        let mut map = PileupChunkHashMap::new("contig1".into());
        let rec = dummy_record(10, Strand::Positive, ModType::FiveMC, 3, 10);
        map.add_record(rec.clone());
        let fetched = map.get_record(10, Strand::Positive, ModType::FiveMC).unwrap();
        assert_eq!(fetched.n_mod, 3);
    }

    #[test]
    fn test_from_chunk() {
        let chunk = PileupChunk {
            reference: "contig1".to_string(),
            records: vec![
                dummy_record(10, Strand::Positive, ModType::FiveMC, 3, 10),
                dummy_record(11, Strand::Negative, ModType::SixMA, 5, 12),
            ],
        };
        let map = PileupChunkHashMap::from_chunk(chunk);
        assert_eq!(map.len(), 2);
        assert!(map.get_record(11, Strand::Negative, ModType::SixMA).is_some());
    }

    #[test]
    fn test_add_records_reference_mismatch_panics() {
        let mut map = PileupChunkHashMap::new("contig1".into());
        let chunk = PileupChunk {
            reference: "wrong_contig".to_string(),
            records: vec![dummy_record(5, Strand::Positive, ModType::FiveMC, 3, 10)],
        };
        let result = map.add_records(chunk);
        assert!(result.is_err());
    }

    #[test]
    fn test_retain_max_mod_per_pos_strand() {
        let mut map = PileupChunkHashMap::new("contig1".into());
        map.add_record(dummy_record(10, Strand::Positive, ModType::FiveMC, 3, 10));
        map.add_record(dummy_record(10, Strand::Positive, ModType::SixMA, 5, 10));
        map.add_record(dummy_record(11, Strand::Positive, ModType::FiveMC, 1, 10));

        map.retain_max_mod_per_pos_strand();

        assert_eq!(map.len(), 2); // 10+ strand should retain only SixMA (5 > 3)
        assert!(map.get_record(10, Strand::Positive, ModType::FiveMC).is_none());
        assert!(map.get_record(10, Strand::Positive, ModType::SixMA).is_some());
        assert!(map.get_record(11, Strand::Positive, ModType::FiveMC).is_some());
    }

    #[test]
    fn test_remove_windows_with_high_methylation() {
        let mut map = PileupChunkHashMap::new("contig1".into());
        for pos in 0..10 {
            map.add_record(dummy_record(pos, Strand::Positive, ModType::FiveMC, 9, 10)); // 0.9 methylation
        }
        map.remove_windows_with_n_high_methylation(5, 0.7, 3);
        assert!(map.len() < 10); // Some records should be removed
    }

    #[test]
    fn test_filter_chunk_coverage_filtering() {
        let mut map = PileupChunkHashMap::new("contig1".into());
        map.add_record(dummy_record(1, Strand::Positive, ModType::FiveMC, 5, 4)); // Low cov
        map.add_record(dummy_record(2, Strand::Positive, ModType::FiveMC, 5, 10)); // High cov
        let _ = filter_chunk(&mut map, 5, 8, 2, 0.4);
        assert!(map.get_record(1, Strand::Positive, ModType::FiveMC).is_none());
        assert!(map.get_record(2, Strand::Positive, ModType::FiveMC).is_some());
    }
}