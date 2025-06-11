// src/cli.rs
use clap::{Parser, ValueEnum};
/// A CLI tool that processes a file with optional numeric parameters.
#[derive(Parser, Debug)]
#[command(name = "mixedmotif", version, about = "Methylation Motif Pairs")]
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
        default_value = "captured_mthylation",
        value_name = "OUT",
        help = "Output file path"
    )]
    pub out: String,


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
