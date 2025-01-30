use std::{collections::HashSet, hash::Hash};

use crate::sequence::{Contig, EqualLengthDNASet, MethylationLevel};
use crate::model::BetaBernoulliModel;
use ahash::{HashMap, HashMapExt};
use anyhow::Result;
use bio::bio_types::annot::contig;
use log::{debug, info};
use ordered_float::OrderedFloat;
use petgraph::{
    algo::min_spanning_tree, graph::{Graph, NodeIndex}, Direction
};
use std::cmp::Ordering;
use std::collections::BinaryHeap;
use utils::{iupac::IupacBase, motif::Motif};
use utils::{modtype::ModType, motif, motif::MotifLike, pileup, strand::Strand};
use itertools::Itertools;


extern crate ordered_float;
pub fn motif_search(
    contig: Contig,
    window_size: usize,
    min_kl_divergence: f64,
    min_base_probability: f64,
    max_branching_with_no_improvement: usize,
) -> Result<(), Box<dyn std::error::Error>> {
    println!("Searching for motifs in contig: {}", contig.reference);

    // Inititate motif result
    let mut motif_result: Vec<motif::Motif> = Vec::new();

    // Iterate over mod types
    for mod_type in ModType::iter() {
        info!("Searching for motifs of type: {:?}", mod_type);
        let seed_motif = motif::Motif::new(
            mod_type.get_iupac_base().to_string(),
            mod_type.to_string(),
            0,
        )?;

        let mut background_contig_seqs = EqualLengthDNASet::new_empty(window_size * 2 + 1)?;
        if let Some(background_motif_indices_fwd) = contig.find_motif_indeces(&seed_motif) {
            let background_contig_seqs_fwd =
                contig.extract_subsequences(&background_motif_indices_fwd, window_size)?;
            let background_contig_set_fwd = EqualLengthDNASet::new(background_contig_seqs_fwd)?;
            background_contig_seqs.add_other(&background_contig_set_fwd)?;
        }
        if let Some(background_motif_indices_rev) =
            contig.find_complement_motif_indeces(&seed_motif)
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

        let mut methylated_motif_seqs = EqualLengthDNASet::new_empty(window_size * 2 + 1)?;
        if let Some(methylated_motif_indices_fwd) =
            contig.methylated_motif_positions(&seed_motif, MethylationLevel::High, Strand::Positive)
        {
            let methylated_motif_seqs_fwd =
                contig.extract_subsequences(&methylated_motif_indices_fwd, window_size)?;
            let methylated_motif_set_fwd = EqualLengthDNASet::new(methylated_motif_seqs_fwd)?;
            methylated_motif_seqs.add_other(&methylated_motif_set_fwd)?;
        }
        if let Some(methylated_motif_indices_rev) =
            contig.methylated_motif_positions(&seed_motif, MethylationLevel::High, Strand::Negative)
        {
            let methylated_motif_seqs_rev =
                contig.extract_subsequences(&methylated_motif_indices_rev, window_size)?;
            let methylated_motif_set_rev =
                EqualLengthDNASet::new(methylated_motif_seqs_rev)?.reverse_complement();
            methylated_motif_seqs.add_other(&methylated_motif_set_rev)?;
        }

        if methylated_motif_seqs.len() == 0 {
            debug!("No methylated positions found for mod type: {:?}", mod_type);
            continue;
        }

        let mut seed_node_methylation_count = MotifMethylationCount::new(seed_motif.clone());
        seed_node_methylation_count.update(&contig);

        let mut keep_searching = true;
        while keep_searching {     
            let mut motif_graph = MotifGraph::new();
            motif_graph.add_node(seed_motif.clone(), 0.0, 0.0);
    
            let mut priority_que = MotifHeap {
                heap: BinaryHeap::new(),
            };
            priority_que.push(motif_graph.get_node(&seed_motif).unwrap().clone());

            let mut best_score = 0.0;
            let mut rounds_since_best_score = 0;

            while let Some(node) = priority_que.pop() {
                rounds_since_best_score += 1;


                let mut node_methylation_count = MotifMethylationCount::new(node.motif.clone());
                node_methylation_count.update(&contig);

                let background_contig_seqs_node =
                    match background_contig_seqs.get_motif_matching_sequences(&node.motif) {
                        Ok(seq) => seq,
                        _ => {
                            debug!(
                                "No matching background sequences found for motif: {:?}",
                                node.motif
                            );
                            continue;
                        }
                    };
                let pssm_background = background_contig_seqs_node.pssm(0.0);

                let methylated_motif_seqs_temp =
                    match methylated_motif_seqs.get_motif_matching_sequences(&node.motif) {
                        Ok(seq) => seq,
                        _ => {
                            debug!(
                                "No matching methylated sequences found for motif: {:?}",
                                node.motif
                            );
                            continue;
                        }
                    };
                let pssm_methylated = methylated_motif_seqs_temp.pssm(0.0);

                let kl_divergence = pssm_methylated.kl_divergence(&pssm_background)?;

                for (i, kl_divergence_value) in kl_divergence.into_iter().enumerate() {
                    if kl_divergence_value > min_kl_divergence {
                        let motif_growth_index = node.motif.position as usize - window_size + i;

                        let valid_bases = pssm_methylated.bases_above_freq(i as u8, min_base_probability);
                        
                        for base in all_unique_combinations(valid_bases)? {
                            let base_str = base.into_iter().map(|b| b.to_string()).collect();
                            let mut new_motif = node.motif.clone();
                            new_motif.grow_motif_at(motif_growth_index as i8, base_str)?;

                            if motif_graph.check_node(&new_motif) {
                                continue;
                            }
                            let mut new_motif_counts = MotifMethylationCount::new(new_motif.clone());
                            new_motif_counts.update(&contig);

                            let score = scoring_function(&new_motif_counts, &node_methylation_count);
                            if score > best_score {
                                best_score = score;
                                rounds_since_best_score = 0;
                            }

                            let priority = original_priority_function(&new_motif_counts, &node_methylation_count, 1.0);
                            if motif_graph.add_node(new_motif.clone(), priority, score) {
                                priority_que.push(motif_graph.get_node(&new_motif).unwrap().clone());
                            } else if motif_graph.get_node(&new_motif).unwrap().score < score {
                                motif_graph.update_node(&new_motif, priority, score);
                            }
                        }
                    }
                }
                
                let best_motif = motif_graph.highest_scoring_node().unwrap();
                info!("Expanded motif: {:<15} score:{:<15.3} | Best: {:<15}, score:{:<15.3}", 
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
            
            if best_score > 0.2 {
                let best_motif = motif_graph.highest_scoring_node().unwrap();
                info!("Keeping Best motif: {:<15} score:{:<15.3}", best_motif.motif.as_string(), best_motif.score);
                motif_result.push(best_motif.motif.clone());
                // Remove this motif from methylated and background sequences
                methylated_motif_seqs = methylated_motif_seqs.remove_motif_matching_sequences(&best_motif.motif)?;
                background_contig_seqs = background_contig_seqs.remove_motif_matching_sequences(&best_motif.motif)?;
            } else {
                keep_searching = false;
            }
        }
        info!("Finished searching for motifs of type: {:?}", mod_type);
        let pretty_motif_result = motif_result.iter().map(|m| m.as_pretty_string()).collect::<Vec<String>>();
        info!("Found {:?} motifs", pretty_motif_result);
    }
    Ok(())
}

struct MotifMethylationCount {
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

    fn update(&mut self, contig: &Contig) {
        self.high_methylated_positions = contig.methylated_motif_count(
            &self.motif,
            MethylationLevel::High,
        );
        self.lowly_methylated_positions = contig.methylated_motif_count(
            &self.motif,
            MethylationLevel::Low,
        );
    }
}

fn scoring_function(motif: &MotifMethylationCount, motif_other: &MotifMethylationCount) -> f64 {
    let beta = BetaBernoulliModel::new_with_params(motif.high_methylated_positions as f64, motif.lowly_methylated_positions as f64);
    let beta_other = BetaBernoulliModel::new_with_params(motif_other.high_methylated_positions as f64, motif_other.lowly_methylated_positions as f64);
    let mean_diff = beta.mean() - beta_other.mean();
    let std_dev = (beta.standard_deviation() + beta_other.standard_deviation()) / 2.0;
    let score = beta.mean() * -std_dev.log10() * mean_diff;
    score


}



fn priority_function(motif: &MotifMethylationCount, motif_other: &MotifMethylationCount) -> f64 {
    // Calculate the priority
    log_odds_ratio(
        motif.high_methylated_positions as f64,
        motif.lowly_methylated_positions as f64,
        motif_other.high_methylated_positions as f64,
        motif_other.lowly_methylated_positions as f64,
        1.0,
    )
}

fn original_priority_function(motif: &MotifMethylationCount, motif_other: &MotifMethylationCount, pseudocount: f64) -> f64 {
    let motif_odds = (motif.high_methylated_positions as f64 + pseudocount) /( motif.lowly_methylated_positions as f64 + pseudocount);
    let motif_other_odds = (motif_other.high_methylated_positions as f64 + pseudocount) / (motif_other.lowly_methylated_positions as f64 + pseudocount);
    let priority = (1.0 - motif_odds) / motif_other_odds;
    priority
}

fn log_odds_ratio(p1: f64, p2: f64, q1: f64, q2: f64, psudo_count: f64) -> f64 {
    let odds_ratop = ((p1 + psudo_count) / (q1 + psudo_count)) / ((p2 + psudo_count) / (q2 + psudo_count));
    -odds_ratop.ln()
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



/// A node in the graph. `score` and `priority` are optional.
#[derive(Debug, Clone, PartialEq)]
pub struct MotifNode {
    pub motif: motif::Motif,
    pub priority: f64,
    pub score: f64,
}

impl MotifNode {
    /// Create a node with *both* priority and score.
    pub fn new(motif: motif::Motif, priority: f64, score: f64) -> Self {
        Self {
            motif,
            priority: priority,
            score: score,
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
    pub fn add_node(&mut self, motif: motif::Motif, priority: f64, score: f64) -> bool {
        if self.node_map.contains_key(&motif) {
            return false;
        }
        let idx = self
            .graph
            .add_node(MotifNode::new(motif.clone(), priority, score));
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

fn all_unique_combinations(items: Vec<IupacBase>) -> Result<Vec<Vec<IupacBase>>> {
    let mut combos = Vec::new();
    for i in 1..=items.len() {
        for combo in items.clone().into_iter().combinations(i) {
            combos.push(combo.clone())
        }
    }
    Ok(combos)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_motif_graph() {
        let mut g = MotifGraph::new();

        let motif_a = motif::Motif::new("GATC", "6mA", 1).unwrap();
        let motif_b = motif::Motif::new("GATC", "6mA", 1).unwrap();
        let motif_c = motif::Motif::new("GATC", "6mA", 1).unwrap();
        let motif_d = motif::Motif::new("GATG", "6mA", 1).unwrap();

        // Add nodes, some with scores, some without
        assert_eq!(g.add_node(motif_b.clone(), 1.0, 2.0), true);
        assert_eq!(g.add_node(motif_b.clone(), 2.0, 3.0), false); // already present
        assert_eq!(g.add_node(motif_c.clone(), 2.0, 3.0), false);
        assert_eq!(g.add_node(motif_d.clone(), 0.0, 0.0), true);

        // Update node that had no score
        assert_eq!(g.update_node(&motif_a, 3.0, 4.0), true);
        // Try updating a non-existent node
        let motif_x = motif::Motif::new("AAAA", "6mA", 1).unwrap();
        assert_eq!(g.update_node(&motif_x, 0.0, 0.0), false);

        // Add edges
        assert_eq!(g.add_edge(&motif_a, &motif_b), true);
        // Duplicate edge
        assert_eq!(g.add_edge(&motif_a, &motif_b), false);
        // Edge from missing node
        assert_eq!(g.add_edge(&motif_x, &motif_b), false);

        // Check neighbors
        let nbrs_a: Vec<&MotifNode> = g.neighbors(&motif_a).unwrap().collect();
        assert_eq!(nbrs_a.len(), 1);
        assert_eq!(nbrs_a[0].motif.mod_type, ModType::SixMA);

        // Highest scoring node
        // - A = score=4.0, B=score=2.0, C=None
        let best = g.highest_scoring_node().unwrap();
        assert_eq!(best.motif.mod_type, ModType::SixMA);
        assert_eq!(best.score, 4.0);
    }
}
