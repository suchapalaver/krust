//! K-mer counting and output.
//!
//! This module provides the main k-mer counting functionality, using parallel
//! processing to efficiently count canonical k-mers across all sequences in a FASTA file.

use crate::{
    cli::OutputFormat,
    kmer::{unpack_to_string, Kmer, KmerLength},
    reader::read,
};
use bytes::Bytes;
use dashmap::DashMap;
use rayon::prelude::{ParallelBridge, ParallelIterator};
use rustc_hash::FxHasher;
use serde::Serialize;
use std::{
    collections::HashMap,
    error::Error,
    fmt::Debug,
    hash::BuildHasherDefault,
    io::{stdout, BufWriter, Error as IoError, Write},
    path::Path,
};
use thiserror::Error;

/// Errors that can occur during k-mer processing.
#[derive(Debug, Error)]
pub enum ProcessError {
    /// Error reading input FASTA file.
    #[error("Unable to read input: {0}")]
    ReadError(#[from] Box<dyn Error>),

    /// Error writing output.
    #[error("Unable to write output: {0}")]
    WriteError(#[from] IoError),

    /// Error serializing JSON output.
    #[error("Unable to serialize JSON: {0}")]
    JsonError(#[from] serde_json::Error),
}

/// A k-mer with its count, used for JSON serialization.
#[derive(Serialize)]
struct KmerCount {
    kmer: String,
    count: i32,
}

/// Counts k-mers in a FASTA file and writes results to stdout.
///
/// Reads sequences from the FASTA file at `path`, counts all canonical k-mers
/// of length `k`, and writes the results to stdout in FASTA-like format.
///
/// # Errors
///
/// Returns `ProcessError::ReadError` if the file cannot be read.
/// Returns `ProcessError::WriteError` if output cannot be written.
pub fn run<P>(path: P, k: usize) -> Result<(), ProcessError>
where
    P: AsRef<Path> + Debug,
{
    run_with_options(path, k, OutputFormat::Fasta, 1)
}

/// Counts k-mers with configurable output format and filtering.
///
/// # Arguments
///
/// * `path` - Path to the FASTA file
/// * `k` - K-mer length
/// * `format` - Output format (fasta, tsv, or json)
/// * `min_count` - Minimum count threshold (k-mers below this are excluded)
///
/// # Errors
///
/// Returns `ProcessError` on read, write, or serialization errors.
pub fn run_with_options<P>(
    path: P,
    k: usize,
    format: OutputFormat,
    min_count: i32,
) -> Result<(), ProcessError>
where
    P: AsRef<Path> + Debug,
{
    let counts = count_kmers(&path, k)?;
    output_counts(counts, format, min_count)
}

/// Counts k-mers and returns them as a HashMap.
///
/// This is the main library API for counting k-mers without writing to stdout.
///
/// # Arguments
///
/// * `path` - Path to the FASTA file
/// * `k` - K-mer length (must be 1-32)
///
/// # Returns
///
/// A HashMap mapping k-mer strings to their counts.
///
/// # Errors
///
/// Returns an error if:
/// - `k` is outside the valid range (1-32)
/// - The file cannot be read
pub fn count_kmers<P>(path: P, k: usize) -> Result<HashMap<String, i32>, Box<dyn Error>>
where
    P: AsRef<Path> + Debug,
{
    // Validate k-mer length upfront to provide a clear error
    let k_len = KmerLength::new(k)?;

    let kmer_map = KmerMap::new().build(read(path)?, k)?;
    Ok(kmer_map.into_hashmap(k_len))
}

fn output_counts(
    counts: HashMap<String, i32>,
    format: OutputFormat,
    min_count: i32,
) -> Result<(), ProcessError> {
    let mut buf = BufWriter::new(stdout());
    let filtered: Vec<_> = counts
        .into_iter()
        .filter(|(_, count)| *count >= min_count)
        .collect();

    match format {
        OutputFormat::Fasta => {
            for (kmer, count) in filtered {
                writeln!(buf, ">{count}\n{kmer}")?;
            }
        }
        OutputFormat::Tsv => {
            for (kmer, count) in filtered {
                writeln!(buf, "{kmer}\t{count}")?;
            }
        }
        OutputFormat::Json => {
            let json_data: Vec<KmerCount> = filtered
                .into_iter()
                .map(|(kmer, count)| KmerCount { kmer, count })
                .collect();
            serde_json::to_writer_pretty(&mut buf, &json_data)?;
            writeln!(buf)?;
        }
    }

    buf.flush()?;
    Ok(())
}

/// A custom `DashMap` w/ `FxHasher`.
type DashFx = DashMap<u64, i32, BuildHasherDefault<FxHasher>>;

struct KmerMap(DashFx);

impl KmerMap {
    fn new() -> Self {
        Self(DashMap::with_hasher(
            BuildHasherDefault::<FxHasher>::default(),
        ))
    }

    fn build(
        self,
        sequences: rayon::vec::IntoIter<Bytes>,
        k: usize,
    ) -> Result<Self, Box<dyn Error>> {
        sequences.for_each(|seq| self.process_sequence(&seq, k));
        Ok(self)
    }

    fn process_sequence(&self, seq: &Bytes, k: usize) {
        if seq.len() < k {
            return;
        }
        let mut i = 0;

        while i <= seq.len() - k {
            let sub = seq.slice(i..i + k);

            match Kmer::from_sub(sub) {
                Ok(unpacked) => {
                    self.process_valid_kmer(unpacked);
                    i += 1;
                }
                Err(err) => {
                    // Skip past the invalid base
                    i += err.position + 1;
                }
            }
        }
    }

    fn process_valid_kmer(&self, unpacked: Kmer) {
        let packed = unpacked.pack();

        // Optimization: check if we've already seen this exact k-mer
        if let Some(mut count) = self.0.get_mut(&packed.packed_bits()) {
            *count += 1;
        } else {
            // Not seen before - compute canonical form
            let canonical = packed.canonical();
            *self.0.entry(canonical.packed_bits()).or_insert(0) += 1;
        }
    }

    fn into_hashmap(self, k: KmerLength) -> HashMap<String, i32> {
        self.0
            .into_iter()
            .par_bridge()
            .map(|(packed_bits, count)| {
                let kmer_string = unpack_to_string(packed_bits, k);
                (kmer_string, count)
            })
            .collect()
    }
}
