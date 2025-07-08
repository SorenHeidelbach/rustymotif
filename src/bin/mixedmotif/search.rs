
use crate::{
    model::{Beta, BetaMixture},
    sequence::{Contig, EqualLengthDNASet, PSSM, MethylationLevel},
};

use ahash::{HashMap, HashMapExt};
use anyhow::Result;
use itertools::Itertools;
use log::{debug, info};
use std::fs::File;
use csv::{WriterBuilder};
use petgraph::{
    graph::{Graph, NodeIndex},
    Direction,
};
use std::{cmp::Ordering, collections::BinaryHeap};
use rustymotif_utils::{
    iupac::IupacBase, modtype::ModType, motif::{self, Motif, MotifLike},
    strand::Strand,
};
use serde::Serialize;

fn all_unique_combinations(bases: Vec<IupacBase>) -> Result<Vec<Vec<IupacBase>>> {
    let mut combinations = Vec::new();
    for i in 1..=bases.len() {
        for combination in bases.iter().cloned().combinations(i) {
            combinations.push(combination);
        }
    }
    Ok(combinations)
}

pub fn motif_search(
    contig: Contig,
    max_branching_with_no_improvement: usize,
    max_low_score_motifs: usize,
    min_score: f64,
    write_intermediate_motifs: Option<&str>,
    window_size: usize,
    min_kl_divergence: f64,
    min_base_probability: f64,
) -> Result<Vec<MotifResult>> {
    println!("Searching for motifs in contig: {}", contig.reference);
    let mut motif_results: Vec<MotifResult> = Vec::new();
    let mut contig_search_clone = contig.clone();

    let _intermediate_result_writer = match write_intermediate_motifs {
        Some(path) => {
            let file = File::create(path)?;
            let mut writer = WriterBuilder::new()
                .delimiter(b'\t')
                .has_headers(false)
                .from_writer(file);
            writer.write_record(
                &[
                    "contig_id",
                    "motif",
                    "motif_seq",
                    "motif_mod_type",
                    "motif_mod_position",
                    "methylated_beta_alpa",
                    "methylated_beta_beta",
                    "methylated_beta_mean",
                    "non_methylated_beta_alpha",
                    "non_methylated_beta_beta",
                    "non_methylated_beta_mean",
                    "methylated_beta_pi",
                    "score",
                ],
            )?;
            Some(writer)
        }
        None => None,
    };

    for mod_type in ModType::iter() {
        info!("Searching for motifs of type: {:?}", mod_type);
        // Initialize seed motif 
        let seed_motif = motif::Motif::new(
            mod_type.get_iupac_base().to_string(),
            mod_type.to_string(),
            0,
        )?;

        // Get records for the mod type
        let mut records = contig_search_clone.get_by_mod_type(mod_type.clone());
        
        // Fit beta to all data to get the background methylation rate
        let mut seed_node_beta = Beta::from_data(&records);
        if seed_node_beta.beta < 1.1 {
            debug!("Seed beta is too low, setting to 1.1");
            seed_node_beta.beta = 1.1;
        }
        let mut seed_model = BetaMixture::new(
            Beta::new(10.0, 1.0),
            seed_node_beta.clone(),
            0.01,
        );
        match seed_model.fit_em_fixed_false(
            &mut records,
            1.0,
            500,
            1e-9
        ) {
            Ok(_) => {
                debug!("Seed model fitted successfully");
            }
            Err(e) => {
                debug!("Failed to fit seed model: {}", e);
                continue;
            }
        }
        
        // Drop the records variable to free the immutable borrow
        drop(records);

        // Create background sequences
        debug!("Creating background sequences");
        let mut background_contig_seqs = EqualLengthDNASet::new_empty(window_size * 2 + 1)?;

        if let Some(background_motif_indices_fwd) = contig.find_motif_indices(&seed_motif)? {
            let background_contig_seqs_fwd =
                contig.extract_subsequences(&background_motif_indices_fwd, window_size)?;
            let background_contig_set_fwd = EqualLengthDNASet::new(background_contig_seqs_fwd)?;
            background_contig_seqs.add_other(&background_contig_set_fwd)?;
        }

        if let Some(background_motif_indices_rev) =
            contig.find_complement_motif_indices(&seed_motif)?
        {
            let background_contig_seqs_rev =
                contig.extract_subsequences(&background_motif_indices_rev, window_size)?;
            let background_contig_seqs_set_rev =
                EqualLengthDNASet::new(background_contig_seqs_rev)?.reverse_complement();
            background_contig_seqs.add_other(&background_contig_seqs_set_rev)?;
        }

        if background_contig_seqs.len() == 0 {
            debug!(
                "No valid positions found for nucleotide: {:?}",
                mod_type.get_iupac_base()
            );
            continue;
        }
        debug!("Created background sequences");

        // Create methylated sequences
        debug!("Creating methylated sequences");
        let mut methylated_motif_seqs = EqualLengthDNASet::new_empty(window_size * 2 + 1)?;

        if let Some(methylated_motif_indices_fwd) =
            contig.methylated_motif_positions(&seed_motif, MethylationLevel::High, Strand::Positive)?
        {
            let methylated_motif_seqs_fwd =
                contig.extract_subsequences(&methylated_motif_indices_fwd, window_size)?;
            if !methylated_motif_seqs_fwd.is_empty() {
                let methylated_motif_set_fwd = EqualLengthDNASet::new(methylated_motif_seqs_fwd)?;
                methylated_motif_seqs.add_other(&methylated_motif_set_fwd)?;
            }
        } else {
            debug!(
                "No methylated positions found on the forward strand for mod type: {:?}",
                mod_type
            );
        }

        let mut keep_searching = true;
        let mut low_score_motifs = 0;
        while keep_searching {
            debug!("Searching for motifs of type: {:?}", mod_type);
            let current_records = contig_search_clone.get_by_mod_type(mod_type.clone());
            debug!("    Found {} records for motif", current_records.len());
            if current_records.len() < 50 {
                break;
            }
            let mut motif_graph = MotifGraph::new();
            // add seed motif
            motif_graph.add_node(
                seed_motif.clone(),
                0.0,
                0.0,
                seed_model.clone(),
                None,
            );

            let mut priority_que = MotifHeap {
                heap: BinaryHeap::new(),
            };
            // push seed to queue
            priority_que.push(motif_graph.get_node(&seed_motif).unwrap().clone());

            let mut best_score = 0.0;
            let mut best_motif = seed_motif.clone();
            let mut best_node = motif_graph.get_node(&seed_motif).unwrap().clone();
            let mut rounds_since_best_score = 0;

            while let Some(node) = priority_que.pop() {
                debug!("Expanding motif: {:25} | Mixing: {:7.4} | Beta mean: {:7.4} | Score: {:7.4}", node.motif.as_pretty_string(), 
                    node.model.pi, 
                    node.model.true_model.mean(),
                    node.score);
                rounds_since_best_score += 1;

                let Ok(background_contig_seqs_node) =
                    background_contig_seqs.get_motif_matching_sequences(&node.motif)
                else {
                    debug!(
                        "No matching background sequences found for motif: {:?}",
                        node.motif
                    );
                    continue;
                };

                let Ok(methylated_motif_seqs_temp) =
                    methylated_motif_seqs.get_motif_matching_sequences(&node.motif)
                else {
                    debug!(
                        "No matching methylated sequences found for motif: {:?}",
                        node.motif
                    );
                    continue;
                };

                let pssm_background = background_contig_seqs_node.pssm(0.0);
                let pssm_methylated = methylated_motif_seqs_temp.pssm(0.0);

                let kl_divergence = pssm_methylated.kl_divergence(&pssm_background)?;

                for i in 0..rounds_since_best_score*2+1 {
                    if 1.0 > min_kl_divergence {
                        let motif_growth_index =
                            node.motif.position as i8 - window_size as i8 + i as i8;

                        let bases_to_expand = vec![IupacBase::A, IupacBase::C, IupacBase::G, IupacBase::T];
                        for base in bases_to_expand {
                            let expansion_base = vec![base.to_string()]; 
                            // Create new motif
                            let mut new_motif = node.motif.clone();
                            new_motif.grow_motif_at(motif_growth_index, expansion_base)?;

                            if motif_graph.has_node(&new_motif) {
                                continue;
                            }

                            // Score the new motif
                            let records = contig_search_clone.get_records_for_motif(&new_motif)?;
                            if records.len() < 10 {
                                continue;
                            }
                            let mut new_model = seed_model.clone();
                            match new_model.fit_em_fixed_false(
                                &records,
                                3.0,
                                100,
                                1e-6,
                            ) {
                                Ok(_) => {}
                                Err(e) => {
                                    debug!("Failed to fit new model: {}", e);
                                    continue;
                                }
                            }

                            let score = new_model.pi * new_model.true_model.mean();
                            let priority = -score;

                            
                            let new_node = MotifNode::new(
                                new_motif.clone(),
                                priority,
                                score,
                                new_model.clone(),
                                Some(pssm_methylated.clone()),
                            );
                            if score > best_score {
                                debug!(
                                    "New best motif: {:?} -> {:?}",
                                    best_node.motif, best_motif
                                );
                                best_score = score;
                                best_motif = new_motif.clone();
                                best_node = new_node.clone();
                                rounds_since_best_score = 0;
                            }
                            priority_que.push(new_node.clone());
                            motif_graph.add_node(
                                new_motif.clone(),
                                priority,
                                score,
                                new_model,
                                Some(pssm_methylated.clone()),
                            );
                            motif_graph.add_edge(&node.motif, &new_motif);
                        }
                    }
                }
                // If no new motifs were found, check if we should stop searching
                if rounds_since_best_score > max_branching_with_no_improvement {
                    debug!(
                        "No new motifs found for {} rounds, stopping search",
                        max_branching_with_no_improvement
                    );
                    break;
                }
            }

            
            // Remove the best node from contig records and methylated positions
            background_contig_seqs.remove_motif_matching_sequences(&best_node.motif)?;
            methylated_motif_seqs.remove_motif_matching_sequences(&best_node.motif)?;
            contig_search_clone.remove_records_for_motif(&best_node.motif)?;
            // Remove records for the best motif from contig

            if best_node.score > min_score {
                info!(
                    "Found motif: {:15} | Score: {:15.4}",
                    best_node.motif.as_pretty_string(),
                    best_node.score
                );
                let motif_result = MotifResult {
                    contig_id: contig.reference.clone(),
                    motif_seq: best_node.motif.sequence_string(),
                    motif_mod_type: mod_type.to_string().to_string(),
                    motif_mod_position: best_node.motif.position as u8,
                    model: best_node.model.clone(),
                };
                motif_results.push(motif_result);
                low_score_motifs = 0;
            } else {
                low_score_motifs += 1;
                if low_score_motifs > max_low_score_motifs {
                    keep_searching = false;
                }
            }
        }
    }
    Ok(motif_results)
}

#[derive(Serialize, Debug, Clone, PartialEq)]
pub struct MotifResult {
    pub contig_id: String,
    pub motif_seq: String,
    pub motif_mod_type: String,
    pub motif_mod_position: u8,
    pub model: BetaMixture,
}

impl MotifResult {
    pub fn new(
        motif: Motif,
        contig: &Contig,
        model: BetaMixture,
    ) -> Result<Self> {
        Ok(Self {
            contig_id: contig.reference.clone(),
            motif_seq: motif.sequence_string(),
            motif_mod_type: motif.mod_type.to_string().to_string(),
            motif_mod_position: motif.position as u8,
            model: model.clone(),
        })
    }
}

/// A node in the graph. `score` and `priority` are optional.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct MotifNode {
    pub motif: motif::Motif,
    pub priority: f64,
    pub score: f64,
    #[serde(flatten)]
    pub model: BetaMixture,
    pub search_pssm: Option<PSSM>,
}

impl MotifNode {
    /// Create a node with *both* priority and score.
    pub fn new(
        motif: motif::Motif,
        priority: f64,
        score: f64,
        model: BetaMixture,
        search_pssm: Option<PSSM>,
    ) -> Self {
        Self {
            motif,
            priority: priority,
            score: score,
            model: model,
            search_pssm,
        }
    }
}

impl Ord for MotifNode {
    fn cmp(&self, other: &Self) -> Ordering {
        other
            .priority
            .partial_cmp(&self.priority)
            .unwrap_or(Ordering::Equal)
    }
}

impl PartialOrd for MotifNode {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Eq for MotifNode {}

struct MotifHeap {
    heap: BinaryHeap<MotifNode>,
}

impl MotifHeap {
    fn push(&mut self, motif: MotifNode) {
        self.heap.push(motif);
    }

    fn pop(&mut self) -> Option<MotifNode> {
        self.heap.pop()
    }
}

pub struct MotifGraph {
    graph: Graph<MotifNode, ()>,
    node_map: HashMap<motif::Motif, NodeIndex>,
}

impl MotifGraph {
    pub fn new() -> Self {
        Self {
            graph: Graph::new(),
            node_map: HashMap::new(),
        }
    }

    // -------- Node Management --------
    pub fn add_node(
        &mut self,
        motif: motif::Motif,
        priority: f64,
        score: f64,
        model: BetaMixture,
        search_pssm: Option<PSSM>,
    ) -> bool {
        if self.node_map.contains_key(&motif) {
            return false;
        }
        let idx = self
            .graph
            .add_node(MotifNode::new(motif.clone(), priority, score, model, search_pssm));
        self.node_map.insert(motif, idx);
        true
    }

    pub fn has_node(&self, motif: &motif::Motif) -> bool {
        self.node_map.contains_key(motif)
    }

    pub fn update_node(&mut self, motif: &motif::Motif, priority: f64, score: f64) -> bool {
        match self.node_map.get(motif) {
            Some(&idx) => {
                let node_data = self.graph.node_weight_mut(idx).expect("NodeIndex invalid");
                node_data.priority = priority;
                node_data.score = score;
                true
            }
            None => false,
        }
    }

    pub fn get_node(&self, motif: &motif::Motif) -> Option<&MotifNode> {
        let &idx = self.node_map.get(motif)?;
        Some(&self.graph[idx])
    }

    // -------- Edge Management --------
    pub fn add_edge(&mut self, from: &motif::Motif, to: &motif::Motif) -> bool {
        let (Some(&idx_from), Some(&idx_to)) = (self.node_map.get(from), self.node_map.get(to))
        else {
            return false;
        };

        if self.graph.find_edge(idx_from, idx_to).is_some() {
            return false;
        }

        self.graph.add_edge(idx_from, idx_to, ());
        true
    }

    // -------- Queries --------
    pub fn neighbors(&self, motif: &motif::Motif) -> Option<impl Iterator<Item = &MotifNode>> {
        let &idx = self.node_map.get(motif)?;
        Some(
            self.graph
                .neighbors_directed(idx, Direction::Outgoing)
                .map(|nbr_idx| &self.graph[nbr_idx]),
        )
    }

    /// Returns the node (by reference) with the highest `score`, ignoring nodes
    /// that have `score=None`. If multiple have the same score, returns the first encountered.
    pub fn highest_scoring_node(&self) -> Option<&MotifNode> {
        self.graph
            .node_indices()
            .map(|idx| {
                let node = &self.graph[idx];
                (node.score, node)
            })
            .max_by(|(score_a, _), (score_b, _)| {
                score_a.partial_cmp(score_b).unwrap_or(Ordering::Equal)
            })
            .map(|(_, node)| node)
    }
}


fn _search_scoring_function(
    _true_beta: &Beta,
    previous_mixture: f64,
    current_mixture: f64,
) -> Result<f64> {
    let score = current_mixture * (current_mixture - previous_mixture);
    Ok(score)
}