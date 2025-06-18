use std::collections::HashMap;
use std::hash::Hash;
use std::path::Path;
use std::fs::File;
use std::io::{BufReader, BufRead};
use anyhow::Result;
use utils::motif::Motif;
use log::{info, debug, warn};

pub fn load_contig_bin(
    contig_bin_path: &Path,
) -> Result<HashMap<String, Vec<String>>> {
    let file = File::open(contig_bin_path)
        .expect("Error opening contig bin file");
    
    let reader = BufReader::new(file);
    let mut bin_contigs = HashMap::new();
    for line in reader.lines() {
        let line = line.expect("Error reading line from contig bin file");
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 2 {
            continue; // Skip lines that do not have at least two parts
        }
        let contig = parts[0].to_string();
        let bin  = parts[1].to_string();
        bin_contigs
            .entry(bin)
            .or_insert_with(Vec::new)
            .push(contig);
    }
    Ok(bin_contigs)


}
pub fn load_bin_motifs(
    bin_motifs_path: &Path,
) -> Result<HashMap<String, Vec<Motif>>> {
    let file = File::open(bin_motifs_path)?;
    
    let reader = BufReader::new(file);
    let mut bin_motifs = HashMap::new();
    for line in reader.lines().skip(1) {
        let line = line?;
        let parts: Vec<&str> = line.split('\t').collect();
        let bin = parts[0].to_string();


        // Forward motif
        let motif_seq = parts[1];
        let motif_position = parts[2].parse::<u8>()?;
        let motif_modtype = parts[3];
        let motif = Motif::new(motif_seq, motif_modtype, motif_position)?;
        bin_motifs
                .entry(bin.clone())
                .or_insert_with(Vec::new)
                .push(motif.clone());
        // Reverse motif
        if parts[7].is_empty() {
            // If there is no reverse motif, we can skip this line
            
            continue;
        }

        let reverse_motif_seq = parts[7];
        let reverse_motif_position = parts[8].parse::<u8>()?;
        let reverse_motif = Motif::new(reverse_motif_seq, motif_modtype, reverse_motif_position)?;

        if reverse_motif.eq(&motif) {
            continue;
        }

        bin_motifs
            .entry(bin)
            .or_insert_with(Vec::new)
            .push(reverse_motif);
    }
    Ok(bin_motifs)
}

