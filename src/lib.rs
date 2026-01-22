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
//! ### Builder API (Recommended)
//!
//! The builder API provides a fluent interface for configuring k-mer counting:
//!
//! ```rust,no_run
//! use kmerust::builder::KmerCounter;
//!
//! fn main() -> Result<(), Box<dyn std::error::Error>> {
//!     // Simple usage
//!     let counts = KmerCounter::new()
//!         .k(21)?
//!         .count("sequences.fa")?;
//!
//!     // With options
//!     let counts = KmerCounter::new()
//!         .k(21)?
//!         .min_count(5)
//!         .count("sequences.fa")?;
//!
//!     for (kmer, count) in counts {
//!         println!("{kmer}: {count}");
//!     }
//!     Ok(())
//! }
//! ```
//!
//! ### Direct API
//!
//! For simpler use cases, the direct API is also available:
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
//!
//! ## Limitations
//!
//! - **K-mer length:** Limited to 1-32 bases (64-bit packing uses 2 bits per base)

#[cfg(feature = "async")]
pub mod async_api;
pub mod builder;
pub mod cli;
pub mod error;
pub mod kmer;
#[cfg(feature = "mmap")]
pub mod mmap;
pub mod progress;
pub(crate) mod reader;
pub mod run;
pub mod streaming;
