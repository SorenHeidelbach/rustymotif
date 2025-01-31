use std::{collections::HashSet, hash::Hash};

use crate::{
    model::BetaBernoulliModel,
    sequence::{Contig, EqualLengthDNASet, MethylationLevel},
};
use ahash::{HashMap, HashMapExt};
use anyhow::Result;
use bio::bio_types::annot::contig;
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
    iupac::IupacBase, modtype::ModType, motif, motif::Motif, motif::MotifLike, pileup,
    strand::Strand,
};
use serde::{Deserialize, Serialize};
extern crate ordered_float;


pub fn motif_search(
    contig: Contig,
    window_size: usize,
    min_kl_divergence: f64,
    min_score: f64,
    min_base_probability: f64,
    max_branching_with_no_improvement: usize,
    max_low_score_motifs: usize,
) -> Result<Vec<MotifResult>> {
    println!("Searching for motifs in contig: {}", contig.reference);

    // Inititate motif result
    let mut motif_results: Vec<MotifResult> = Vec::new();

    // Iterate over mod types
    for mod_type in ModType::iter() {
        info!("Searching for motifs of type: {:?}", mod_type);
        let seed_motif = motif::Motif::new(
            mod_type.get_iupac_base().to_string(),
            mod_type.to_string(),
            0,
        )?;

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

        debug!("Creating methylated sequences");
        let mut methylated_motif_seqs = EqualLengthDNASet::new_empty(window_size * 2 + 1)?;

        if let Some(methylated_motif_indices_fwd) =
            contig.methylated_motif_positions(&seed_motif, MethylationLevel::High, Strand::Positive)?
        {
            let methylated_motif_seqs_fwd =
                contig.extract_subsequences(&methylated_motif_indices_fwd, window_size)?;
            let methylated_motif_set_fwd = EqualLengthDNASet::new(methylated_motif_seqs_fwd)?;
            methylated_motif_seqs.add_other(&methylated_motif_set_fwd)?;
        } else {
            debug!(
                "No methylated positions found on the forward strand for mod type: {:?}",
                mod_type
            );
        }

        if let Some(methylated_motif_indices_rev) =
            contig.methylated_motif_positions(&seed_motif, MethylationLevel::High, Strand::Negative)?
        {
            let methylated_motif_seqs_rev =
                contig.extract_subsequences(&methylated_motif_indices_rev, window_size)?;
            let methylated_motif_set_rev =
                EqualLengthDNASet::new(methylated_motif_seqs_rev)?.reverse_complement();
            methylated_motif_seqs.add_other(&methylated_motif_set_rev)?;
        } else {
            debug!(
                "No methylated positions found on the reverse strand for mod type: {:?}",
                mod_type
            );
        }

        if methylated_motif_seqs.len() == 0 {
            debug!("No methylated positions found for mod type: {:?}", mod_type);
            continue;
        }
        debug!("Created methylated sequences");

        debug!("Initializing seed node");
        let mut seed_node_methylation_count = MotifMethylationCount::new(seed_motif.clone());
        seed_node_methylation_count.update(&contig)?;


        debug!("Starting search");
        let mut keep_searching = true;
        let priority_calculator = PriorityCalculator::new(1.0);
        let score_calculator = ScoreCalculator::new(1.0);

        let mut low_score_motifs = 0;
        while keep_searching {
            let mut motif_graph = MotifGraph::new();
            motif_graph.add_node(
                seed_motif.clone(),
                0.0,
                0.0,
                seed_node_methylation_count.clone(),
            );

            let mut priority_que = MotifHeap {
                heap: BinaryHeap::new(),
            };
            priority_que.push(motif_graph.get_node(&seed_motif).unwrap().clone());

            let mut best_score = 0.0;
            let mut rounds_since_best_score = 0;

            while let Some(node) = priority_que.pop() {
                rounds_since_best_score += 1;

                let mut node_methylation_count = MotifMethylationCount::new(node.motif.clone());
                node_methylation_count.update(&contig)?;

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

                for (i, kl_divergence_value) in kl_divergence.into_iter().enumerate() {
                    if kl_divergence_value > min_kl_divergence {
                        let motif_growth_index = node.motif.position as usize - window_size + i;

                        let valid_bases =
                            pssm_methylated.bases_above_freq(i as u8, min_base_probability);

                        for base in all_unique_combinations(valid_bases)? {
                            // Create new motif
                            let base_str = base.into_iter().map(|b| b.to_string()).collect();
                            let mut new_motif = node.motif.clone();
                            new_motif.grow_motif_at(motif_growth_index as i8, base_str)?;

                            let mut new_motif_counts =
                                MotifMethylationCount::new(new_motif.clone());
                            new_motif_counts.update(&contig)?;

                            let score = score_calculator
                                .ad_hoc(&new_motif_counts, &node_methylation_count)?;
                            let priority = priority_calculator
                                .priority_ad_hoc(&new_motif_counts, &seed_node_methylation_count);
                            if score > best_score {
                                best_score = score;
                                rounds_since_best_score = 0;
                            }

                            if motif_graph.add_node(
                                new_motif.clone(),
                                priority,
                                score,
                                new_motif_counts,
                            ) {
                                priority_que
                                    .push(motif_graph.get_node(&new_motif).unwrap().clone());
                            } else if motif_graph.get_node(&new_motif).unwrap().score < score {
                                motif_graph.update_node(&new_motif, priority, score);
                            }
                            motif_graph.add_edge(&node.motif, &new_motif);
                        }
                    }
                }

                let best_motif = motif_graph.highest_scoring_node().unwrap();
                info!(
                    "Expanded motif: {:<15} score:{:<15.3} | Best: {:<15}, score:{:<15.3}",
                    node.motif.as_pretty_string(),
                    node.score,
                    best_motif.motif.as_pretty_string(),
                    best_motif.score
                );

                if rounds_since_best_score > max_branching_with_no_improvement {
                    debug!("No improvement in {max_branching_with_no_improvement} rounds, stopping search");
                    break;
                }
            }

            if best_score > min_score {
                let best_motif = motif_graph.highest_scoring_node().unwrap();
                info!(
                    "Keeping Best motif: {:<15} score:{:<15.3}",
                    best_motif.motif.as_string(),
                    best_motif.score
                );


                // Remove this motif from methylated and background sequences
                methylated_motif_seqs =
                    methylated_motif_seqs.remove_motif_matching_sequences(&best_motif.motif)?;
                background_contig_seqs =
                    background_contig_seqs.remove_motif_matching_sequences(&best_motif.motif)?;
                
                let motif_result = MotifResult::new(
                    best_motif.motif.clone(),
                    &contig,
                    best_motif.score.clone(),
                )?;
                motif_results.push(motif_result);
            } else {
                low_score_motifs += 1;
                if low_score_motifs > max_low_score_motifs {
                    debug!("Too many low score motifs, stopping search");
                    keep_searching = false;
                }
            }
        }
        info!("Finished searching for motifs of type: {:?}", mod_type);
        let pretty_motif_results = motif_results
            .iter()
            .map(|m| m.motif.as_pretty_string())
            .collect::<Vec<String>>();
        info!("Found {:?} motifs", pretty_motif_results);
    }
    Ok(motif_results)
}

#[derive(Serialize, Debug, Clone, PartialEq)]
pub struct MotifResult {
    contig_id: String,
    motif: Motif,
    motif_seq: String,
    motif_mod_type: String,
    motif_mod_position: u8,
    low_count: usize,
    middle_count: usize,
    high_count: usize,
    score: f64,
}

impl MotifResult {
    pub fn new(
        motif: Motif,
        contig: &Contig,
        score: f64,
    ) -> Result<Self> {
        let low_count = contig.methylated_motif_count(&motif, MethylationLevel::Low)?;
        let middle_count = contig.methylated_motif_count(&motif, MethylationLevel::Middle)?;
        let high_count = contig.methylated_motif_count(&motif, MethylationLevel::High)?;
        Ok(Self {
            contig_id: contig.reference.clone(),
            motif: motif.clone(),
            motif_seq: motif.sequence_string(),
            motif_mod_type: motif.mod_type.to_string().to_string(),
            motif_mod_position: motif.position.clone(),
            low_count,
            middle_count,
            high_count,
            score,
        })
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct MotifMethylationCount {
    motif: motif::Motif,
    high_methylated_positions: usize,
    lowly_methylated_positions: usize,
}

impl MotifMethylationCount {
    fn new(motif: motif::Motif) -> Self {
        Self {
            motif,
            high_methylated_positions: 0,
            lowly_methylated_positions: 0,
        }
    }

    fn update(&mut self, contig: &Contig) -> Result<()> {
        self.high_methylated_positions =
            contig.methylated_motif_count(&self.motif, MethylationLevel::High)?;
        self.lowly_methylated_positions =
            contig.methylated_motif_count(&self.motif, MethylationLevel::Low)?;
        Ok(())
    }
}

/// A node in the graph. `score` and `priority` are optional.
#[derive(Debug, Clone, PartialEq)]
pub struct MotifNode {
    pub motif: motif::Motif,
    pub priority: f64,
    pub score: f64,
    pub motif_info: MotifMethylationCount,
}

impl MotifNode {
    /// Create a node with *both* priority and score.
    pub fn new(
        motif: motif::Motif,
        priority: f64,
        score: f64,
        count: MotifMethylationCount,
    ) -> Self {
        Self {
            motif,
            priority: priority,
            score: score,
            motif_info: count,
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
        count: MotifMethylationCount,
    ) -> bool {
        if self.node_map.contains_key(&motif) {
            return false;
        }
        let idx = self
            .graph
            .add_node(MotifNode::new(motif.clone(), priority, score, count));
        self.node_map.insert(motif, idx);
        true
    }

    pub fn check_node(&self, motif: &motif::Motif) -> bool {
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

#[derive(Debug, Clone, PartialEq)]
struct ScoreCalculator {
    pseudocount: f64,
}

impl ScoreCalculator {
    fn new(pseudocount: f64) -> Self {
        Self { pseudocount }
    }

    fn score_log_odds_ratio(
        &self,
        motif: &MotifMethylationCount,
        motif_other: &MotifMethylationCount,
    ) -> Result<f64> {
        let score = -log_odds_ratio(
            motif.high_methylated_positions as f64,
            motif.lowly_methylated_positions as f64,
            motif_other.high_methylated_positions as f64,
            motif_other.lowly_methylated_positions as f64,
            self.pseudocount,
        );
        Ok(score)
    }

    fn score_log_odds_ratio_z_score(
        &self,
        motif: &MotifMethylationCount,
        motif_other: &MotifMethylationCount,
    ) -> Result<f64> {
        let score = -log_odds_ratio_z_score(
            motif.high_methylated_positions as f64,
            motif.lowly_methylated_positions as f64,
            motif_other.high_methylated_positions as f64,
            motif_other.lowly_methylated_positions as f64,
            self.pseudocount,
        );
        Ok(score)
    }

    fn ad_hoc(
        &self,
        motif: &MotifMethylationCount,
        motif_other: &MotifMethylationCount,
    ) -> Result<f64> {
        let mut beta = BetaBernoulliModel::new_with_params(
            motif.high_methylated_positions as f64,
            motif.lowly_methylated_positions as f64,
        );
        beta.update(self.pseudocount as usize, self.pseudocount as usize);
        let mut beta_other = BetaBernoulliModel::new_with_params(
            motif_other.high_methylated_positions as f64,
            motif_other.lowly_methylated_positions as f64,
        );
        beta_other.update(self.pseudocount as usize, self.pseudocount as usize);

        let mean_diff = beta.mean() - beta_other.mean();
        let std_dev = (beta.standard_deviation() + beta_other.standard_deviation()) / 2.0;
        let score = beta.mean() * -std_dev.log10() * mean_diff;
        Ok(score)
    }
}

#[derive(Debug, Clone, PartialEq)]
struct PriorityCalculator {
    pseudocount: f64,
}

impl PriorityCalculator {
    fn new(pseudocount: f64) -> Self {
        Self { pseudocount }
    }

    fn priority_ad_hoc(
        &self,
        motif: &MotifMethylationCount,
        motif_other: &MotifMethylationCount,
    ) -> f64 {
        let motif_odds = (motif.high_methylated_positions as f64 + self.pseudocount)
            / (motif.lowly_methylated_positions as f64 + self.pseudocount);
        let motif_other_odds = (motif_other.high_methylated_positions as f64 + self.pseudocount)
            / (motif_other.lowly_methylated_positions as f64 + self.pseudocount);
        let priority = (1.0 - motif_odds) / motif_other_odds;
        priority
    }

    fn priority_log_odds_ratio(
        &self,
        motif: &MotifMethylationCount,
        motif_other: &MotifMethylationCount,
    ) -> f64 {
        log_odds_ratio(
            motif.high_methylated_positions as f64,
            motif.lowly_methylated_positions as f64,
            motif_other.high_methylated_positions as f64,
            motif_other.lowly_methylated_positions as f64,
            self.pseudocount,
        )
    }

    fn priority_log_odds_ratio_z_score(
        &self,
        motif: &MotifMethylationCount,
        motif_other: &MotifMethylationCount,
    ) -> f64 {
        log_odds_ratio_z_score(
            motif.high_methylated_positions as f64,
            motif.lowly_methylated_positions as f64,
            motif_other.high_methylated_positions as f64,
            motif_other.lowly_methylated_positions as f64,
            self.pseudocount,
        )
    }
}

fn log_odds_ratio(p1: f64, p2: f64, q1: f64, q2: f64, psudo_count: f64) -> f64 {
    let odds_ratio =
        ((p1 + psudo_count) / (q1 + psudo_count)) / ((p2 + psudo_count) / (q2 + psudo_count));
    -odds_ratio.ln()
}

fn log_odds_ratio_z_score(p1: f64, p2: f64, q1: f64, q2: f64, psudo_count: f64) -> f64 {
    let p1_pseudo = p1 + psudo_count;
    let p2_pseudo = p2 + psudo_count;
    let q1_pseudo = q1 + psudo_count;
    let q2_pseudo = q2 + psudo_count;

    let log_odds_ratio = log_odds_ratio(p1, p2, q1, q2, psudo_count);
    // handle inf case
    if log_odds_ratio.is_infinite() {
        return log_odds_ratio;
    }
    let var = (1.0 / p1_pseudo) + (1.0 / p2_pseudo) + (1.0 / q1_pseudo) + (1.0 / q2_pseudo);
    let z_score = log_odds_ratio / var.sqrt();
    z_score
}

fn all_unique_combinations(items: Vec<IupacBase>) -> Result<Vec<Vec<IupacBase>>> {
    let mut combos = Vec::new();
    for i in 1..=items.len() {
        for combo in items.clone().into_iter().combinations(i) {
            combos.push(combo.clone())
        }
    }
    Ok(combos)
}
