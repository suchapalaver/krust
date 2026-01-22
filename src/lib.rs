//! # kmerust
//!
//! A fast, parallel [k-mer](https://en.wikipedia.org/wiki/K-mer) counter for DNA sequences in FASTA files.
//!
//! ## Features
//!
//! - Parallel processing using [rayon](https://docs.rs/rayon) and [dashmap](https://docs.rs/dashmap)
//! - Outputs canonical k-mers (lexicographically smaller of k-mer and reverse complement)
//! - Supports k-mer lengths from 1 to 32
//! - Handles sequences with N bases (skips invalid k-mers)
//! - Compatible output format with [Jellyfish](https://github.com/gmarcais/Jellyfish)
//!
//! ## CLI Usage
//!
//! ```bash
//! # Count 21-mers in a FASTA file
//! kmerust 21 sequences.fa > kmers.txt
//!
//! # Count 5-mers
//! kmerust 5 sequences.fa > kmers.txt
//! ```
//!
//! ## Output Format
//!
//! Output is written to stdout in FASTA-like format:
//! ```text
//! >{count}
//! {canonical_kmer}
//! ```
//!
//! ## Library Usage
//!
//! ```rust,no_run
//! use kmerust::run::count_kmers;
//! use std::path::PathBuf;
//!
//! fn main() -> Result<(), Box<dyn std::error::Error>> {
//!     let path = PathBuf::from("sequences.fa");
//!     let counts = count_kmers(&path, 21)?;
//!     for (kmer, count) in counts {
//!         println!("{kmer}: {count}");
//!     }
//!     Ok(())
//! }
//! ```

pub mod cli;
pub mod error;
pub mod kmer;
#[cfg(feature = "mmap")]
pub mod mmap;
pub mod progress;
pub(crate) mod reader;
pub mod run;
pub mod streaming;
