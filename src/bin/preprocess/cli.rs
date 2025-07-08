// src/cli.rs
use clap::{Parser, ValueEnum};
/// A CLI tool that processes a file with optional numeric parameters.
#[derive(Parser, Debug)]
#[command(name = "Preprocess", version, about = "Preprocesses methyltion pileup data for downstream usage")]
pub struct Cli {
    #[arg(
        value_name = "PILEUP",
        help = "File path to the pileup file with methylation data"
    )]
    pub pileup: String,

    #[arg(
        long,
        short,
        default_value = "preprocess_test3.bed",
        value_name = "OUT",
        help = "Output file path"
    )]
    pub out: String,

    
    
    #[arg(
        long,
        default_value = "5",
        help = "Minimum coverage required to consider a position"
    )]
    pub min_cov: u32,
    #[arg(long, default_value = "4", help = "Window size to search for motifs")]
    pub window_size: usize,

    #[arg(
        long,
        default_value = "0.4",
        help = "Fraction threshold for high methylation density"
    )]
    pub fraction_threshold: f32,

    #[arg(
        long,
        default_value = "2",
        help = "Number of high methylation positions required in a window"
    )]
    pub n_high_methylation: usize,

    #[arg(long, short, default_value = "5", help = "Number of threads to use")]
    pub threads: u32,


    #[arg(
        value_enum,
        long,
        default_value = "normal",
        value_name = "VERBOSITY",
        help = "Verbosity level"
    )]
    pub verbosity: LogLevel,

}

#[derive(Debug, clap::ValueEnum, Clone, Copy, PartialEq, Eq)]
pub enum LogLevel {
    Verbose,
    Normal,
    Silent,
}

impl std::fmt::Display for LogLevel {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            LogLevel::Verbose => write!(f, "verbose"),
            LogLevel::Normal => write!(f, "normal"),
            LogLevel::Silent => write!(f, "silent"),
        }
    }
}
