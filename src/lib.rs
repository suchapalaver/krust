//! # krust
//!
//! `krust` is a [k-mer](https://en.wikipedia.org/wiki/K-mer) counter
//! written in Rust and run from the command line that will output canonical
//! k-mers and their frequency across the records in a fasta file.
//!
//! `krust` prints to `stdout`, writing, on alternate lines:  
//! ```>{frequency}```  
//! ```{canonical k-mer}```  
//!
//! `krust` uses [`dashmap`](https://docs.rs/crate/dashmap/4.0.2),
//! [`rust-bio`](https://docs.rs/bio/0.38.0/bio/), [`rayon`](https://docs.rs/rayon/1.5.1/rayon/), 
//! and [`fxhash`](https://crates.io/crates/fxhash).
//!
//! Run `krust` on the test data in the [`krust` Github repo](https://github.com/suchapalaver/krust),
//! searching for kmers of length 5, like this:  
//! ```$ cargo run --release 5 cerevisae.pan.fa > output.tsv```  
//! or, searching for kmers of length 21:  
//! ```$ cargo run --release 21 cerevisae.pan.fa > output.tsv```  
//!
//! Future:
//! /// Returns k-mer counts for individual sequences in a fasta file.  
//! ```fn single_sequence_canonical_kmers(filepath: String, k: usize) {}```      

use bio::{alphabets::dna::revcomp, io::fasta};
use dashmap::DashMap;
use fxhash::FxHasher;
use rayon::prelude::*;
use std::{cmp::min, env, hash::BuildHasherDefault};

/// This struct needs to be public as part of the way `krust` works,
/// a simple struct for parsing the command line arguments we need.
pub struct Config {
    pub kmer_len: usize,
    pub filepath: String,
}

impl Config {
    pub fn new(mut args: env::Args) -> Result<Config, &'static str> {
        args.next();

        let kmer_len = match args.next() {
            Some(arg) => arg.parse().unwrap(),
            None => return Err("Problem with k-mer length input"),
        };
        let filepath = args.next().unwrap();

        Ok(Config { kmer_len, filepath })
    }
}

/// A custom `DashMap` w/ `FxHasher`.
pub type KrustMap = DashMap<Box<[u8]>, u64, BuildHasherDefault<FxHasher>>;

/// Reads sequences from fasta records in parallel using [`rayon`](https://docs.rs/rayon/1.5.1/rayon/).
/// Using [`Dashmap`](https://docs.rs/dashmap/4.0.2/dashmap/struct.DashMap.html) allows updating a
/// single hashmap in parallel.
/// Ignores substrings containing `N`.
/// Canonicalizes by lexicographically smaller of k-mer/reverse-complement.
/// Returns a hashmap of canonical k-mers (keys) and their frequency in the data (values).
pub fn canonicalize_kmers(filepath: String, k: usize) -> Result<KrustMap, &'static str> {
    let canonical_hash: KrustMap = DashMap::with_hasher(BuildHasherDefault::<FxHasher>::default());

    fasta::Reader::from_file(&filepath)
        .unwrap()
        .records()
        .into_iter()
        .par_bridge()
        .for_each(|record| {
            let seq: &[u8] = record.as_ref().unwrap().seq();

            for i in 0..(seq.len() + 1).saturating_sub(k) {
                if !&seq[i..i + k].contains(&b'N') {
                    let substring: &[u8] = &seq[i..i + k];
                    let canonical_kmer = Box::from(min(substring, &revcomp(substring)));

                    *canonical_hash.entry(canonical_kmer).or_insert(0) += 1;
                } else {
                }
            }
        });

    Ok(canonical_hash)
}
