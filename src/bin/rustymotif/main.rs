use clap::Parser;
use env_logger::Env;
use log::info;
use std::path::Path;
use rustymotif_utils::motif;
use rustymotif_utils::pileup;
use anyhow::Result;

mod cli;
mod data;
mod fasta_reader;
mod model;
mod motif_discovery;
mod search;
mod sequence;
fn main() -> Result<()> {
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
    motif_discovery::rustymotif(&args)?;
    info!("Finished motif discovery");
    Ok(())
}
