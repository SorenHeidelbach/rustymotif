use clap::builder::styling::AnsiColor;
use clap::Parser;
use env_logger::Env;
use log::{info, debug, warn};
use std::clone;
use std::collections::{HashMap, HashSet};
use std::path::Path;
use utils::motif;
use utils::pileup;
use utils::modtype::ModType;
use utils::strand::Strand;
use std::fs::File;
use std::time::Instant;
use csv::{WriterBuilder};
use anyhow::Result;
use crate::data::AnnotatedPileupRecord;
use crate::io::load_bin_motifs;
use crate::io::load_contig_bin;

mod cli;
mod data;
mod fasta_reader;
mod sequence;
mod io;
fn main() {
    let args = cli::Cli::parse();
    // Set up logging level
    match args.verbosity {
        cli::LogLevel::silent => {
            env_logger::Builder::from_env(Env::default().default_filter_or("off")).init();
        }
        cli::LogLevel::normal => {
            env_logger::Builder::from_env(Env::default().default_filter_or("info")).init();
        }
        cli::LogLevel::verbose => {
            env_logger::Builder::from_env(Env::default().default_filter_or("debug")).init();
        }
    }

    // Create output directory
    info!("Running motif methylation state");
    let out_path = Path::new(&args.out);
    match out_path.exists() {
        true => {
            panic!("Output directory already exists");
        }
        false => match std::fs::create_dir(out_path) {
            Ok(_) => info!("Created output directory"),
            Err(e) => panic!("Could not create output directory: {}", e),
        },
    }

    // Run the main function
    captured_methylation(&args)
        .unwrap_or_else(|e| {
            eprintln!("Error: {}", e);
            std::process::exit(1);
        });
}



fn captured_methylation(args: &cli::Cli) -> Result<()> {
    let global_timer = Instant::now();

    info!("Starting motif methylation state processing");
    info!("Loading reference file: {}", args.reference);
    let reference_file = Path::new(&args.reference);
    let reference = fasta_reader::read_fasta_file(reference_file)
        .map_err(|e| anyhow::anyhow!("Error reading reference file: {}", e))?;
    info!("Loaded {} reference records", reference.len());

    let bed_file = File::open(&args.bed)
        .map_err(|e| anyhow::anyhow!("Could not open bed file: {} ({})", args.bed, e))?;
    let mut pileup_reader = pileup::PileupChunkReader::new(bed_file, 0);
   
    let contig_bin_path = Path::new(&args.contig_bin);
    info!("Loading contig bin from: {}", args.contig_bin);
    let bin_contig_map = load_contig_bin(contig_bin_path)?;

    let bin_motifs_path = Path::new(&args.bin_motifs);
    info!("Loading bin motifs from: {}", args.bin_motifs);
    let mut bin_motifs = load_bin_motifs(bin_motifs_path)
        .map_err(|e| anyhow::anyhow!("Error loading bin motifs: {}", e))?;
    debug!("Loaded {} bin motifs", bin_motifs.values().len());


    // Create a map to hold motifs for each contig
    let mut contig_motif_map: HashMap<String, Vec<motif::Motif>> = HashMap::new();
    for (bin, motifs) in bin_motifs.iter_mut() {
        if let Some(contigs) = bin_contig_map.get(bin) {
            for contig in contigs {
                contig_motif_map.entry(contig.clone())
                    .or_insert_with(Vec::new)
                    .extend(motifs.iter().cloned());
            }
        } else {
            warn!("No contigs found for bin: {}", bin);
        }
    }

    let outdir = Path::new(&args.out);
    let output = outdir.join(format!("motif_annotated_pileup.bed"));
    let output_string = output.to_str()
        .ok_or_else(|| anyhow::anyhow!("Could not convert path to string"))?;
    
    let result_file = File::create(output_string)?;
    let mut result_writer = WriterBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_writer(result_file);
    result_writer.write_record(&[
        "contig_id",
        "mod_type",
        "position",
        "strand",
        "n_mod",
        "n_canonical",
        "n_valid_cov",
        "n_diff",                    
        "in_motif",
        "motifs",
    ])?;

    let summary_output = outdir.join(format!("motif_annotated_pileup_summary.tsv"));
                        let summary_output_string = summary_output.to_str()
                            .ok_or_else(|| anyhow::anyhow!("Could not convert path to string"))?;
                        let mut summary_writer = WriterBuilder::new()
                            .delimiter(b'\t')
                            .has_headers(true)
                            .from_path(summary_output_string)?;
                        summary_writer.write_record(&[
                            "contig_id",
                            "mod_type",
                            "n_motifs",
                            "motifs",
                            "n_records",
                            "n_in_motif",
                            "n_not_in_motif",
                            "n_mod_in_motif",
                            "n_canonical_in_motif",
                            "n_valid_cov_in_motif",
                            "n_diff_in_motif",
                            "n_gt_threshold_in_motif",
                            "n_mod_not_in_motif",
                            "n_canonical_not_in_motif",
                            "n_valid_cov_not_in_motif",
                            "n_diff_not_in_motif",
                            "n_gt_threshold_not_in_motif",
                        ])?;
    
    info!("Processing pileup file: {}", args.bed);
    loop {
        info!("Processing a batch");
        let timer = Instant::now();
        let chunks = pileup_reader.load_n_chunks(1);
        match chunks {
            Some(chunks) => {
                info!("Loaded batch {:?}", timer.elapsed());
                let mut builder = data::GenomeWorkSpaceBuilder::new();
                for chunk in chunks {
                    let contig_id = &chunk.reference;
                    info!("Processing contig: {}", contig_id);
                    debug!("Adding contig to workspace");
                    let contig = sequence::Contig::new(
                        contig_id,
                        reference.get(contig_id).unwrap().as_str(),
                    );
                    builder.add_contig(contig);
                    debug!("Adding records to contig");
                    builder.push_records(chunk);
                }
                let genome_work_space = builder.build();

                
                
                for ( contig_id, contig) in genome_work_space.contigs.iter() {
                    debug!("Contig ID: {}", contig_id);
                    let motifs = contig_motif_map.get(contig_id)
                        .cloned()
                        .unwrap_or_else(|| vec![]);
                    let motifs_set: HashSet<motif::Motif> = HashSet::from_iter(motifs.iter().cloned());

                    // Iterate over pileup records in the contig
                    for mod_type in ModType::iter() {
                        debug!("Processing mod_type: {}", mod_type);
                        let records = contig.get_by_mod_type(mod_type.clone());
                        if records.is_empty() {
                            continue;
                        }

                        let mut annotated_records = records.iter()
                            .map(|record| {
                                AnnotatedPileupRecord::new(
                                    record.clone().clone(),
                                    false, // in_motif will be set later
                                    vec![], // Clone motifs for each record
                                )
                            })
                            .collect::<Vec<AnnotatedPileupRecord>>();
                        
        
                        // Add motif to records where motif occours
                        for motif in motifs_set.iter() {
                            if motif.mod_type != mod_type {
                                continue;
                            }
                            debug!("Processing motif: {}", motif.as_pretty_string());
                            let output = outdir.join(format!("{}_motif_{}.tsv", contig_id, motif.as_pretty_string()));
                            let output_string = output.to_str()
                            .ok_or_else(|| anyhow::anyhow!("Could not convert path to string"))?;
                        
                            // Get index of where motif occours
                            let motif_index = contig.find_motif_indices(motif)?;
                            let motif_index_set = match motif_index {
                                None => None,
                                Some(indices) => Some(indices.iter().cloned().collect::<HashSet<usize>>()),
                            };

                            let motif_index_reverse = contig.find_complement_motif_indices(motif)?;
                            let motif_index_set_reverse = match motif_index_reverse {
                                None => None,
                                Some(indices) => Some(indices.iter().cloned().collect::<HashSet<usize>>()),
                            };


                            for record in annotated_records.iter_mut() {
                                match record.pileup_record.strand {
                                    Strand::Positive => {
                                        match &motif_index_set {
                                            None => continue,
                                            Some(set) => {
                                                if set.contains(&record.pileup_record.position) {
                                                    record.in_motif = true;
                                                    record.motifs.push(motif.clone());
                                                }
                                            }
                                        }
                                    },
                                    Strand::Negative => {
                                        match &motif_index_set {
                                            None => continue,
                                            Some(set) => {
                                                if set.contains(&record.pileup_record.position) {
                                                    record.in_motif = true;
                                                    record.motifs.push(motif.clone());
                                                }
                                            }
                                        }
                                    }
                                }
                            }
    
                        }
                        // Write annotated records to file
                        for record in &annotated_records {
                            result_writer.write_record(&[
                                &contig_id,
                                &record.pileup_record.mod_type.to_string().to_string(),
                                &record.pileup_record.position.to_string(),
                                &record.pileup_record.strand.to_string(),
                                &record.pileup_record.n_mod.to_string(),
                                &record.pileup_record.n_canonical.to_string(),
                                &record.pileup_record.n_valid_cov.to_string(),
                                &record.pileup_record.n_diff.to_string(),
                                &record.in_motif.to_string(),
                                &record.motifs.iter().map(|m| m.as_pretty_string()).collect::<Vec<String>>().join(","),
                            ])?;
                        }

                        // Summarise values over the contig and write to file
                        let filtered_annotated_records = annotated_records.iter()
                            .filter(|r| r.pileup_record.n_valid_cov > args.summary_min_cov)
                            .collect::<Vec<&AnnotatedPileupRecord>>();

                        let filtered_annotated_records_in_motif = filtered_annotated_records.iter()
                            .filter(|r| r.in_motif)
                            .cloned()
                            .collect::<Vec<&AnnotatedPileupRecord>>();
                        let filtered_annotated_records_not_in_motif = filtered_annotated_records.iter()
                            .filter(|r| !r.in_motif)
                            .cloned()
                            .collect::<Vec<&AnnotatedPileupRecord>>();
                        
                        let n_motifs = motifs_set.len();
                        let n_records = annotated_records.len();
                        let n_in_motif = filtered_annotated_records_in_motif.iter()
                            .count();
                        let n_not_in_motif = n_records - n_in_motif;
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
                            .filter(|r| (r.pileup_record.n_mod / r.pileup_record.n_valid_cov) as f32 > args.summary_min_methylation)
                            .count();
                        let n_mod_not_in_motif = annotated_records.iter()
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
                            .filter(|r| (r.pileup_record.n_mod / r.pileup_record.n_valid_cov) as f32 > args.summary_min_methylation)
                            .count();
                        summary_writer.write_record(&[
                            contig_id,
                            &mod_type.to_string().to_string(),
                            &n_motifs.to_string(),
                            &motifs_set.iter()
                                .filter(|m| m.mod_type == mod_type)
                                .map(|m| m.as_pretty_string())
                                .collect::<Vec<String>>().join(","),
                            &n_records.to_string(),
                            &n_in_motif.to_string(),
                            &n_not_in_motif.to_string(),
                            &n_mod_in_motif.to_string(),
                            &n_canonical_in_motif.to_string(),
                            &n_valid_cov_in_motif.to_string(),
                            &n_diff_in_motif.to_string(),
                            &n_gt_threshold_in_motif.to_string(),
                            &n_mod_not_in_motif.to_string(),
                            &n_canonical_not_in_motif.to_string(),
                            &n_valid_cov_not_in_motif.to_string(),
                            &n_diff_not_in_motif.to_string(),
                            &n_gt_threshold_not_in_motif.to_string(),
                        ])?;
                        info!("Wrote summary for contig: {} mod_type: {}", contig_id, mod_type);
                        info!("Wrote {} records to file: {}", annotated_records.len(), output_string);
                    }   
                }
            }
            None => {
                info!("Contig did not contain any records");
            }
        }

        info!("Finished batch in {:?}", timer.elapsed());

        if pileup_reader.eof_reached {
            break;
        }
    }
    info!("Finished processing in {:?}", global_timer.elapsed());
    Ok(())
}