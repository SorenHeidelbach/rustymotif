// src/cli.rs
use clap::{Parser, ValueEnum};
/// A CLI tool that processes a file with optional numeric parameters.
#[derive(Parser, Debug)]
#[command(name = "rustymotif", version, about = "Methylation Motif Pairs")]
pub struct Cli {
    #[arg(
        value_name = "REFERENCE",
        help = "File path to the fasta file with references"
    )]
    pub reference: String,

    #[arg(
        value_name = "PILEUP",
        help = "File path to the pileup file with methylation data"
    )]
    pub pileup: String,

    #[arg(
        long,
        short,
        default_value = "rustymotif",
        value_name = "OUT",
        help = "Output file path"
    )]
    pub out: String,

    #[arg(
        long,
        default_value = "0.1",
        help = "Middle threshold for methylation"
    )]
    pub middle_threshold: f64,

    #[arg(
        long,
        default_value = "0.75",
        help = "High threshold for methylation"
    )]
    pub high_threshold: f64,

    #[arg(long, default_value = "10", help = "Window size to search for motifs")]
    pub window_size: usize,

    #[arg(
        long,
        default_value = "0.1",
        help = "Minimum KL divergence to consider a motif"
    )]
    pub min_kl_divergence: f64,

    #[arg(
        long,
        default_value = "5",
        help = "Minimum coverage required to consider a position"
    )]
    pub min_cov: u32,

    #[arg(long, short, default_value = "5", help = "Number of threads to use")]
    pub threads: u32,

    #[arg(
        long,
        default_value = "100",
        help = "Number of contigs to load and process at once"
    )]
    pub batch_size: u32,

    #[arg(
        value_enum,
        long,
        default_value = "normal",
        value_name = "VERBOSITY",
        help = "Verbosity level"
    )]
    pub verbosity: LogLevel,
}

#[derive(ValueEnum, Clone, Debug)]
pub enum LogLevel {
    verbose,
    normal,
    silent,
}
