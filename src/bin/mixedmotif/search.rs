use std::{collections::HashSet, hash::Hash};

use crate::{
    model::{
        fit_true_beta_model_em, Beta, BetaBernoulliModel, BetaMixture
    },
    sequence::{Contig, EqualLengthDNASet},
};
use ahash::{HashMap, HashMapExt};
use anyhow::Result;
use bio::stats::bayesian::model;
use strum::IntoEnumIterator;
use strum_macros::EnumIter;
use itertools::{concat, Itertools};
use log::{debug, info};
use ordered_float::OrderedFloat;

use std::{fs::File, path::Path, time::Instant};
use csv::{WriterBuilder};
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


pub fn motif_search(
    contig: Contig,
    max_branching_with_no_improvement: usize,
    max_low_score_motifs: usize,
    min_score: f64,
    write_intermediate_motifs: Option<&str>,
) -> Result<Vec<MotifResult>> {
    println!("Searching for motifs in contig: {}", contig.reference);
    let mut motif_results: Vec<MotifResult> = Vec::new();

    let mut intermediate_result_writer = match write_intermediate_motifs {
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
                    "non_methylated_beta_alpha",
                    "non_methylated_beta_beta",
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
        let seed_motif = motif::Motif::new(
            mod_type.get_iupac_base().to_string(),
            mod_type.to_string(),
            0,
        )?;
        let mut records = contig.get_by_mod_type(mod_type.clone());

        let seed_node_beta = Beta::from_data(&records);
        let mut seed_model = BetaMixture::new(
            Beta::new(5.0, 2.0),
            seed_node_beta.clone(),
            0.05,
        );
        seed_model.fit_em_fixed_false(
            &records,
            1.0,
            100,
            1e-6
        )?;

        let mut keep_searching = true;
        let mut low_score_motifs = 0;
        while keep_searching {
            let mut motif_graph = MotifGraph::new();
            // add seed motif
            motif_graph.add_node(
                seed_motif.clone(),
                0.0,
                0.0,
                seed_model.clone()
            );

            let mut priority_que = MotifHeap {
                heap: BinaryHeap::new(),
            };
            // push seed to queue
            priority_que.push(motif_graph.get_node(&seed_motif).unwrap().clone());

            let mut best_score = 0.0;
            let mut best_motif = seed_motif.clone();
            let mut rounds_since_best_score = 0;

            while let Some(node) = priority_que.pop() {
                rounds_since_best_score += 1;
                let mut expanded_motifs = Vec::new();
                // Expand motif
                for iupac in IupacBase::iter().take(4) {
                    let mut expanded_motif = node.motif.clone();
                    expanded_motif.add_base_upstream(iupac)?;
                    expanded_motifs.push(expanded_motif.clone());

                    let mut expanded_motif = node.motif.clone();
                    expanded_motif.add_base_downstream(iupac)?;
                    expanded_motifs.push(expanded_motif.clone());
                }
                for expanded_motif in expanded_motifs {
                    if motif_graph.check_node(&expanded_motif) {
                        debug!("Expanded motif already in graph, skipping");
                        continue;
                    }
                    // Expanded motif locations
                    let expanded_motif_locations = match contig.find_motif_indices(&expanded_motif.clone())? {
                        Some(locations) => {
                            if locations.is_empty() {
                                continue;
                            }
                            locations
                        }
                        None => {
                            continue;
                        }
                        
                    };
                    let expanded_motif_locations_set: HashSet<_> = expanded_motif_locations
                        .iter().cloned().collect();
                    
                    

                    // Get the records for the expanded motif
                    let expanded_motif_records = &records
                        .iter()
                        .filter(|record| {
                            expanded_motif_locations_set.contains(&record.position)
                        })
                        .cloned()
                        .collect::<Vec<_>>();

                    if expanded_motif_records.len() < 10 {
                        debug!("Too few records found for expanded motif");
                        continue;
                    }

                    // fit new beta model
                    let mut expanded_motif_beta_mixture = BetaMixture::new(
                        Beta::new(5.0, 2.0),
                        seed_node_beta.clone(),
                        0.05,
                    );
                    expanded_motif_beta_mixture.fit_em_fixed_false(
                        &expanded_motif_records,
                        1.0,
                        100,
                        1e-6
                    )?;

                    

                    let expanded_motif_score = search_scoring_function(
                        &expanded_motif_beta_mixture.true_model, 
                        node.model.pi, 
                        expanded_motif_beta_mixture.pi)?;
                    motif_graph.add_node(
                        expanded_motif.clone(), 
                        -expanded_motif_score, 
                        expanded_motif_score, 
                        expanded_motif_beta_mixture.clone());
                    priority_que.push(
                        motif_graph.get_node(&expanded_motif).unwrap().clone()
                    );

                    // Write expanded motif to file


                    if let Some(ref mut writer) = intermediate_result_writer {
                        let record = (
                            contig.reference.clone(),
                            expanded_motif.clone(),
                            expanded_motif.sequence_string(),
                            expanded_motif.mod_type.to_string(),
                            expanded_motif.position,
                            expanded_motif_beta_mixture.true_model.alpha,
                            expanded_motif_beta_mixture.true_model.beta,
                            expanded_motif_beta_mixture.false_model.alpha,
                            expanded_motif_beta_mixture.false_model.beta,
                            expanded_motif_beta_mixture.pi,
                            expanded_motif_score,
                        );
                        writer.serialize(record)?;
                        writer.flush()?;
                    }

                    if expanded_motif_score > best_score {
                        best_score = expanded_motif_score;
                        best_motif = expanded_motif.clone();
                        debug!("New best motif: {:?}", best_motif);
                        rounds_since_best_score = 0;
                    }
                }



                // Grow current motif

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

            let best_motif = motif_graph.highest_scoring_node().unwrap();
            // filter away records at motif indices
            debug!("Removing motif locations from records");
            let best_motif_locations = match contig.find_motif_indices(&best_motif.motif.clone())? {
                Some(locations) => {
                    if locations.is_empty() {
                        continue;
                    }
                    locations
                }
                None => {
                    continue;
                }
                
            };
            let best_motif_locations_set: HashSet<_> = best_motif_locations
                .iter().cloned().collect();
            records.retain(|record| {
                !best_motif_locations_set.contains(&record.position)
            });
            
            if best_score > min_score {
                info!(
                    "Keeping Best motif: {:<15} score:{:<15.3}",
                    best_motif.motif.as_string(),
                    best_motif.score
                );


                let motif_result = MotifResult::new(
                    best_motif.motif.clone(),
                    &contig,
                    best_motif.model.clone(),
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
    model: BetaMixture,
}

impl MotifResult {
    pub fn new(
        motif: Motif,
        contig: &Contig,
        model: BetaMixture,
    ) -> Result<Self> {
        Ok(Self {
            contig_id: contig.reference.clone(),
            motif: motif.clone(),
            motif_seq: motif.sequence_string(),
            motif_mod_type: motif.mod_type.to_string().to_string(),
            motif_mod_position: motif.position.clone(),
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
}

impl MotifNode {
    /// Create a node with *both* priority and score.
    pub fn new(
        motif: motif::Motif,
        priority: f64,
        score: f64,
        model: BetaMixture,
    ) -> Self {
        Self {
            motif,
            priority: priority,
            score: score,
            model: model,
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
    ) -> bool {
        if self.node_map.contains_key(&motif) {
            return false;
        }
        let idx = self
            .graph
            .add_node(MotifNode::new(motif.clone(), priority, score, model));
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


fn search_scoring_function(
    true_beta: &Beta,
    previous_mixture: f64,
    current_mixture: f64,
) -> Result<f64> {
    let score = current_mixture * (current_mixture - previous_mixture);
    Ok(score)
}