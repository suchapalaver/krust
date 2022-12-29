//! # krust
//!
//! `krust` is a [k-mer](https://en.wikipedia.org/wiki/K-mer) counter
//! written in Rust and run from the command line that will output canonical
//! k-mers and their frequency across the records in a fasta file.
//!
//! `krust` prints to `stdout`, writing, on alternate lines, just like the popular [`Jellyfish`](https://github.com/gmarcais/Jellyfish) k-mer counter:
//! ```>{frequency}```
//! ```{canonical k-mer}```
//!
//! `krust` has been tested throughout production against [`jellyfish`](https://github.com/gmarcais/Jellyfish)'s results for the same data sets.
//!
//! `krust` uses [`dashmap`](https://docs.rs/crate/dashmap/4.0.2),
//! [`rust-bio`](https://docs.rs/bio/0.38.0/bio/), [`rayon`](https://docs.rs/rayon/1.5.1/rayon/),
//! and [`fxhash`](https://crates.io/crates/fxhash).
//!
//! Run `krust` on the test data in the [`krust` Github repo](https://github.com/suchapalaver/krust),
//! searching for kmers of length 5, like this:
//! ```$ cargo run --release 5 path/to/cerevisae.pan.fa > output.tsv```
//! or, searching for kmers of length 21:
//! ```$ cargo run --release 21 path/to/cerevisae.pan.fa > output.tsv```
//!
//! Future:
//! - ```fn single_sequence_canonical_kmers(filepath: String, k: usize) {}```
//! Returns k-mer counts for individual sequences in a fasta file.
//! - Testing!

pub mod config;
pub(crate) mod kmer;
pub(crate) mod reader;
pub mod startup;
