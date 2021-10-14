//! # krust
//!
//! `krust` is a [k-mer](https://en.wikipedia.org/wiki/K-mer) counter written in Rust and run from the command line that will output canonical k-mers and their frequency across the records in a fasta file.
//!
//! `krust` prints to `stdout`, writing, on alternate lines:  
//! ```>{frequency}```  
//! ```{canonical k-mer}```  
//!
//! `krust` uses [`rust-bio`](https://docs.rs/bio/0.38.0/bio/), [`rayon`](https://docs.rs/rayon/1.5.1/rayon/), and [`dashmap`](https://docs.rs/crate/dashmap/4.0.2).  
//!
//! Run `krust` on the test data in the [`krust` Github repo](https://github.com/suchapalaver/krust), searching for kmers of length 5, like this:  
//! ```$ cargo run --release 5 cerevisae.pan.fa > output.tsv```  
//! or, searching for kmers of length 21:  
//! ```$ cargo run --release 21 cerevisae.pan.fa > output.tsv```  
//!
//! Future:  
//! A function like fn single_sequence_canonical_kmers(filepath: String, k: usize) {}    
//! Would returns k-mer counts for individual sequences in a fasta file.     

use bio::{alphabets::dna::revcomp, io::fasta};
use dashmap::DashMap;
use rayon::prelude::*;
use std::{cmp::min, env};

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

/// Reads sequences from fasta records in parallel using [`rayon`](https://docs.rs/rayon/1.5.1/rayon/).
/// Using [`Dashmap`](https://docs.rs/dashmap/4.0.2/dashmap/struct.DashMap.html) allows updating single hashmap in parallel.
/// Ignores substrings containing `N`.
/// Canonicalizes by lexicographically smaller of k-mer/reverse-complement.
/// Returns a hashmap of canonical k-mers (keys) and their frequency in the data (values).
pub fn canonicalize_kmer(
    filepath: String,
    k: usize,
) -> Result<DashMap<Box<[u8]>, u64>, &'static str> {
    let canonical_hash: DashMap<Box<[u8]>, u64> = DashMap::new();

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

                    *canonical_hash
                        .entry(Box::from(min(substring, &revcomp(substring))))
                        .or_insert(0) += 1;
                } else {
                }
            }
        });

    Ok(canonical_hash)
}
