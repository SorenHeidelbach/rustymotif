use clap::builder::styling::AnsiColor;
use clap::Parser;
use env_logger::Env;
use log::{info, debug};
use std::clone;
use std::collections::HashMap;
use std::path::Path;
use utils::motif;
use utils::pileup;
use utils::modtype::ModType;
use std::fs::File;
use std::time::Instant;
use csv::{WriterBuilder};
use crate::data::AnnotatedPileupRecord;
use polar::prelude::*;

mod cli;
mod data;
mod fasta_reader;
mod sequence;

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
    match captured_methylation(&args) {
        Ok(_) => info!("Finished running motif methylation state"),
        Err(e) => panic!("Error running motif methylation state: {}", e),
    }
}



fn captured_methylation(args: &cli::Cli) -> Result<(), Box<dyn std::error::Error>> {
    let global_timer = Instant::now();
    info!("Starting motif methylation state processing");
    let reference_file = Path::new(&args.reference);
    let reference = fasta_reader::read_fasta_file(reference_file)
        .map_err(|e| anyhow::anyhow!("Error reading reference file: {}", e))?;
    info!("Loaded {} reference records", reference.len());

    let bed_file = File::open(&args.bed)
        .map_err(|e| anyhow::anyhow!("Could not open bed file: {} ({})", args.bed, e))?;
    let mut pileup_reader = pileup::PileupChunkReader::new(bed_file, 0);

    let bin_motifs = Path::new(&args.bin_motifs);
    if !bin_motifs.exists() {
        return Err(anyhow::anyhow!("Bin motifs file does not exist: {}", args.bin_motifs));
    }
    info!("Loading bin motifs from: {}", args.bin_motifs);
    let bin_motifs = LazyFrame::scan_csv(
        bin_motifs,
        CsvReadOptions::new()
            .has_header(true)
            .with_delimiter(b'\t')
            .with_rechunk(true),
    );
    print!(bin_motifs);

    
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

                
                // Write motif_beta_true, motif_mixture, motif_log_likelihood, motif_log_likelihood_false to a file
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
                for ( contig_id, contig) in genome_work_space.contigs.iter() {
                    debug!("Contig ID: {}", contig_id);

                    let motifs = vec![
                        motif::Motif::new("RGATCY", "4mC", 4)?,
                        motif::Motif::new("RGATCY", "4mC", 4)?,
                        motif::Motif::new("GATC", "6mA", 1)?,
                    ]; 

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
                        for motif in motifs.iter() {
                            if motif.mod_type != mod_type {
                                continue;
                            }
                            debug!("Processing motif: {}", motif.as_pretty_string());
                            let output = outdir.join(format!("{}_motif_{}.tsv", contig_id, motif.as_pretty_string()));
                            let output_string = output.to_str()
                                .ok_or_else(|| anyhow::anyhow!("Could not convert path to string"))?;

                            // Get index of where motif occours
                            let motif_index = contig.find_motif_indices(motif)?;
                            if motif_index.is_none() {
                                debug!("Motif ({}) not found in contig ({})", motif.as_pretty_string(), contig_id);
                                continue;
                            }



                            // Add to pileup records where motif occurs 
                            let motif_index = motif_index.unwrap();
                            let motif_index_set: std::collections::HashSet<usize> = motif_index.iter().cloned().collect();
                            for record in annotated_records.iter_mut() {
                                if motif_index_set.contains(&record.pileup_record.position) {
                                    record.in_motif = true;
                                    record.motifs.push(motif.clone());
                                }
                            }
    
                        }
                        // Write annotated records to file
                        for record in annotated_records {
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