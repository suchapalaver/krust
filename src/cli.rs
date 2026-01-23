//! Command-line interface definition.

use clap::{Parser, ValueEnum};
use std::path::PathBuf;

use crate::input::Input;

/// A fast, parallel k-mer counter for DNA sequences in FASTA files.
///
/// Supports reading from files or stdin (use `-` or omit path for stdin).
///
/// # Examples
///
/// ```bash
/// # Count k-mers from a file
/// krust 21 sequences.fa
///
/// # Count k-mers from stdin
/// cat sequences.fa | krust 21
/// cat sequences.fa | krust 21 -
///
/// # With gzip
/// zcat large.fa.gz | krust 21 > counts.tsv
/// ```
#[derive(Parser, Debug)]
#[command(name = "kmerust")]
#[command(version, author, about, long_about = None)]
pub struct Args {
    /// K-mer length (1-32)
    #[arg(value_parser = parse_k)]
    pub k: usize,

    /// Input file path (use '-' or omit for stdin)
    #[arg(default_value = "-")]
    pub path: PathBuf,

    /// Output format
    #[arg(short, long, value_enum, default_value = "fasta")]
    pub format: OutputFormat,

    /// Minimum count threshold (k-mers below this are excluded)
    #[arg(short, long, default_value = "1")]
    pub min_count: u64,

    /// Suppress informational output (only output k-mer counts)
    #[arg(short, long)]
    pub quiet: bool,
}

impl Args {
    /// Returns the input source (file or stdin).
    #[must_use]
    pub fn input(&self) -> Input {
        Input::from_path(&self.path)
    }
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
