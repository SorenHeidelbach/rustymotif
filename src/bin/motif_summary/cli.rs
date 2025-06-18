// src/cli.rs
use clap::{Parser, ValueEnum};
#[derive(Parser, Debug)]
#[command(name = "motif_summary", version, about = "Summarise a pileup file over nanomotif motifs")]
pub struct Cli {
    #[arg(
        long,
        short = 'r',
        value_name = "REFERENCE",
        help = "File path to the fasta file with references"
    )]
    pub reference: String,

    #[arg(
        long,
        short='p',
        value_name = "METHYLBED",
        help = "File path to the pileup file with methylation data"
    )]
    pub bed: String,

    #[arg(
        long,
        short='m',
        value_name = "BIN_MOTIFS",
        help = "File path to the bin motifs file"
    )]
    pub bin_motifs: String,

    #[arg(
        long,
        short='c',
        value_name = "CONTIG_BIN",
        help = "File path to the contig bin file"
    )]
    pub contig_bin: String,

    #[arg(
        long,
        short='o',
        default_value = "motif_summary",
        value_name = "OUT",
        help = "Output path"
    )]
    pub out: String,

    #[arg(
        long,
        value_name = "summary_min_cov",
        default_value = "5",
        help = "Minimum coverage for summary counts"
    )]
    pub summary_min_cov: u32,

    #[arg(
        long,
        value_name = "summary_min_methylation",
        default_value = "0.7",
        help = "Minimum methylation level for summary counts"
    )]
    pub summary_min_methylation: f32,

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
