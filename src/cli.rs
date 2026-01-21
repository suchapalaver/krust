//! Command-line interface definition.

use clap::{Parser, ValueEnum};
use std::path::PathBuf;

/// A fast, parallel k-mer counter for DNA sequences in FASTA files.
#[derive(Parser, Debug)]
#[command(name = "krust")]
#[command(version, author, about, long_about = None)]
pub struct Args {
    /// K-mer length (1-32)
    #[arg(value_parser = parse_k)]
    pub k: usize,

    /// Path to a FASTA file
    pub path: PathBuf,

    /// Output format
    #[arg(short, long, value_enum, default_value = "fasta")]
    pub format: OutputFormat,

    /// Minimum count threshold (k-mers below this are excluded)
    #[arg(short, long, default_value = "1")]
    pub min_count: i32,

    /// Suppress informational output (only output k-mer counts)
    #[arg(short, long)]
    pub quiet: bool,
}

/// Output format for k-mer counts.
#[derive(Debug, Clone, Copy, ValueEnum, Default)]
pub enum OutputFormat {
    /// FASTA-like format (>{count}\n{kmer})
    #[default]
    Fasta,
    /// Tab-separated values (kmer\tcount)
    Tsv,
    /// JSON array format
    Json,
}

fn parse_k(s: &str) -> Result<usize, String> {
    let k: usize = s
        .parse()
        .map_err(|_| format!("'{s}' is not a valid number"))?;
    if k == 0 {
        return Err("k-mer length must be at least 1".to_string());
    }
    if k > 32 {
        return Err("k-mer length must be at most 32".to_string());
    }
    Ok(k)
}
