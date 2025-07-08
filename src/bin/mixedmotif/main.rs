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
mod score_motifs;
use crate::cli::LogLevel;
#[derive(Parser)]
#[clap(author, version, about, long_about = None)]
struct Args {
    #[clap(short, long, default_value = "output")]
    out: String,
    #[clap(short, long, default_value_t = LogLevel::Normal)]
    verbosity: LogLevel,
    #[clap(long, default_value = "100")]
    max_branching: usize,
    #[clap(long, default_value = "100")]
    max_low_score_motifs: usize,
    #[clap(long, default_value = "25.0")]
    min_score: f64,
    #[clap(long)]
    write_intermediate_motifs: Option<String>,
    #[clap(long, default_value = "10")]
    window_size: usize,
    #[clap(long, default_value = "0.2")]
    min_kl_divergence: f64,
    #[clap(long, default_value = "0.1")]
    min_base_probability: f64,
}

fn main() -> Result<()> {
    let args = cli::Cli::parse();
    // Set up logging level
    match args.verbosity {
        cli::LogLevel::Silent => {
            env_logger::Builder::from_env(Env::default().default_filter_or("off")).init();
        }
        cli::LogLevel::Normal => {
            env_logger::Builder::from_env(Env::default().default_filter_or("info")).init();
        }
        cli::LogLevel::Verbose => {
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
        false => std::fs::create_dir(out_path)?,
    }

    // Run the main function
    motif_discovery::rustymotif(&args)?;
    Ok(())
}
