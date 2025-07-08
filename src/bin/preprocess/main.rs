use rustymotif_utils::preprocess;
use clap::Parser;
use rustymotif_utils::pileup;
use env_logger::Env;
use log::{info, debug, warn};
use rustymotif_utils::pileup::PileupRecord;
use io::{create_methyl_bed_writer, write_records_iter};
use cli::Cli;
use anyhow::Result;

mod cli;
mod io;
fn main() {
    env_logger::Builder::from_env(Env::default().default_filter_or("debug")).init();
    let args = cli::Cli::parse();


    print!("Running preprocess...\n");
    debug!("Running preprocess with arguments: {:?}", args);
    preprocess(args).unwrap_or_else(|e| {
        eprintln!("Error during preprocessing: {}", e);
        std::process::exit(1);
    });
    info!("Finished processing all contigs.");
}



fn preprocess(args: Cli) -> Result<()>  {
    // Initiaizing output writer
    let path = std::path::Path::new(&args.out);
    let mut record_writer = create_methyl_bed_writer(path, args.bgzip)?;
    
    // Initializing pileup reader
    let pileup_reader = pileup::PileupTabixReader::new(
        args.pileup.clone(),
        args.min_cov,
        args.threads as usize
    )?;
    let contig_ids = pileup_reader.reader.seqnames();


    for contig_id in &contig_ids {
        let time = std::time::Instant::now();
        let mut chunk = preprocess::load_contig(
            contig_id, 
            args.pileup.clone(),
            args.threads as usize
        ).unwrap_or_else(|e| {
            warn!("Error loading contig {}, {}\nskipping...", contig_id, e);
            rustymotif_utils::preprocess::PileupChunkHashMap::new(contig_id.to_string())
        });
        info!("Processing contig: {}", contig_id);
        if preprocess::filter_chunk(
            &mut chunk, 
            args.min_cov,
            args.window_size,
            args.n_high_methylation,
            args.fraction_threshold,
        ).is_err() {
            warn!("Skipping contig {} due to processing error", contig_id);
            continue;
        }
        let mut values: Vec<PileupRecord> = chunk.records.into_values().collect();
        values.sort_unstable_by(|r1, r2| r1.position.cmp(&r2.position));
        if write_records_iter(&mut *record_writer, &values).is_err() {
            warn!("Error writing records for contig {}", contig_id);
            continue;
        }
        info!("Wrote records to file for contig {}", contig_id);
        let time_taken = time.elapsed();
        println!("Time taken to read all records: {:?}", time_taken);
    }
    Ok(())
}