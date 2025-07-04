use std::os::unix::raw::time_t;
use clap::Parser;
use utils::motif;
use utils::pileup;
use env_logger::Env;
use log::{info, debug, warn};
use crate::cli::Cli;

mod data;
mod sequence;
mod cli;
fn main() {
    env_logger::Builder::from_env(Env::default().default_filter_or("debug")).init();

    let args = cli::Cli::parse();

    print!("Running preprocess...\n");
    let mut pileup_reader = pileup::PileupTabixReader::new(args.pileup, 5, 20).unwrap();
    let contig_ids = pileup_reader.reader.seqnames();
    for contig_id in &contig_ids {
        println!("Contig ID: {}", contig_id);
        println!("Fetching contig: {}", contig_id);
        let time = std::time::Instant::now();
        pileup_reader.fetch_contig(contig_id).unwrap();
        println!("Reading all records from contig: {}", contig_id);
        match pileup_reader.read_all_from_contig() {
            Some(chunk) => {
                println!("Read {} records from contig {}", chunk.records.len(), contig_id);
            },
            None => {
                println!("No records from contig {}", contig_id);
            }
        }
        let time_taken = time.elapsed();
        println!("Time taken to read all records: {:?}", time_taken);
    }
}

