//! Command-line interface definition.

use clap::{Parser, Subcommand, ValueEnum};
use std::path::PathBuf;

use crate::format::SequenceFormat;
use crate::input::Input;

/// A fast, parallel k-mer counter for DNA sequences in FASTA and FASTQ files.
///
/// Supports reading from files or stdin (use `-` or omit path for stdin).
/// Input format is auto-detected from file extension, or can be specified explicitly.
///
/// # Examples
///
/// ```bash
/// # Count k-mers from a FASTA file
/// krust 21 sequences.fa
///
/// # Count k-mers from a FASTQ file
/// krust 21 reads.fq
///
/// # Count k-mers from stdin
/// cat sequences.fa | krust 21
/// cat sequences.fa | krust 21 -
///
/// # Stdin with explicit format
/// cat reads.fq | krust 21 --input-format fastq
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

    /// Input file format (auto-detected from extension if not specified)
    #[arg(short = 'i', long = "input-format", value_enum, default_value = "auto")]
    pub input_format: SequenceFormat,

    /// Save k-mer counts to index file for later querying
    #[arg(long)]
    pub save: Option<PathBuf>,

    /// Minimum base quality score (Phred, 0-93) for FASTQ filtering.
    /// K-mers containing bases below this threshold are skipped.
    /// Ignored for FASTA input.
    #[arg(short = 'Q', long = "min-quality")]
    pub min_quality: Option<u8>,
}

impl Args {
    /// Returns the input source (file or stdin).
    #[must_use]
    pub fn input(&self) -> Input {
        Input::from_path(&self.path)
    }

    /// Returns the resolved input format.
    ///
    /// If `input_format` is `Auto`, detects from the file extension.
    /// For stdin without explicit format, defaults to FASTA.
    #[must_use]
    pub fn resolved_input_format(&self) -> SequenceFormat {
        self.input_format.resolve(Some(&self.path))
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
    /// Histogram format (count\tfrequency) - count of counts
    Histogram,
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

/// Top-level CLI structure supporting both counting and querying.
#[derive(Parser, Debug)]
#[command(name = "kmerust")]
#[command(
    version,
    author,
    about = "A fast, parallel k-mer counter for DNA sequences"
)]
pub struct Cli {
    #[command(subcommand)]
    pub command: Option<Command>,
}

/// Available subcommands.
#[derive(Subcommand, Debug)]
pub enum Command {
    /// Query k-mer counts from a pre-built index
    Query(QueryArgs),
}

/// Arguments for the query command.
#[derive(Parser, Debug)]
pub struct QueryArgs {
    /// Path to the k-mer index file (.kmix)
    pub index: PathBuf,

    /// K-mer sequence to query (e.g., ACGTACGT)
    pub kmer: String,
}
