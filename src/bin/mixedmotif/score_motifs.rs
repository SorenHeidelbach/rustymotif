use std::{collections::HashSet, fs::File, hash::Hash};

use crate::{
    model::{
        fit_true_beta_model_em, Beta, BetaBernoulliModel, BetaMixture
    },
    sequence::{Contig, EqualLengthDNASet},
};
use ahash::{HashMap, HashMapExt};
use anyhow::{Ok, Result};
use strum::IntoEnumIterator;
use strum_macros::EnumIter;
use itertools::Itertools;
use log::{debug, info};
use ordered_float::OrderedFloat;
use petgraph::{
    algo::min_spanning_tree,
    graph::{Graph, NodeIndex},
    Direction,
};
use std::{cmp::Ordering, collections::BinaryHeap};
use utils::{
    iupac::{self, IupacBase}, modtype::ModType, motif::{self, Motif, MotifLike}, pileup,
    strand::Strand,
}; 
use serde::{Deserialize, Serialize};
extern crate ordered_float;


pub fn score_motif(
    motif: Motif,
    contig: Contig,
    write_motif_records: Option<String>,
) -> Result<Option<BetaMixture>> {
    let mod_type = motif.mod_type;
    debug!("Scoring motif: {:?}", motif);
    let records = contig.get_by_mod_type(mod_type.clone());
    let seed_node_beta = Beta::from_data(&records);

    // Expanded motif locations
    debug!("    Finding motif locations");
    let motif_index = match contig.find_motif_indices(&motif.clone())?  {
        Some(motif_index) => motif_index,
        None => {
            debug!("Motif ({}) not found in contig ({})", motif.as_pretty_string(), contig.reference);
            return Ok(None);
        }
    };
    let motif_index_set: HashSet<usize> = motif_index.iter().cloned().collect();
    let motif_reverse_index = match contig.find_complement_motif_indices(&motif)? {
        Some(motif_reverse_index) => motif_reverse_index,
        None => {
            debug!("Motif ({}) not found in reverse for contig ({})", motif.as_pretty_string(), contig.reference);
            return Ok(None);
        }
    };
    let motif_reverse_index_set: HashSet<usize> = motif_reverse_index.iter().cloned().collect();
    
    debug!("Getting records for motif");
    let motif_records = records
        .iter()
        .filter(|record| {
            motif_index_set.contains(&record.position) && record.strand == utils::strand::Strand::Positive
        })
        .cloned()
        .collect::<Vec<_>>();

    if motif_records.is_empty() {
        debug!("No records found for motif");
        return Ok(None);
    }
    debug!("    Found {} records for motif", motif_records.len());
    if let Some(write_motif_records) = write_motif_records {
        debug!("    Writing motif records to file");
        let motif_records_file = File::create(write_motif_records)?;
        let mut motif_records_writer = csv::WriterBuilder::new()
            .delimiter(b'\t')
            .has_headers(true)
            .from_writer(motif_records_file);
        for record in motif_records.iter() {
        motif_records_writer.serialize(record)?;
        }
        motif_records_writer.flush()?;
    }

    let mut beta_mixture_model = BetaMixture::new(
        Beta::new(5.0, 2.0),
        seed_node_beta.clone(),
        0.05,
    );
    beta_mixture_model.fit_em_fixed_false(&motif_records, 2.0, 100, 1e-6)?;
    debug!("Mean of motif Beta: {:?}", beta_mixture_model.true_model.mean());
    debug!("Standard deviuation of motyif Beta: {:?}", beta_mixture_model.true_model.standard_deviation());

    Ok(Some(beta_mixture_model))
}
        
        



  

