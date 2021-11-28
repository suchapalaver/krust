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
//! ```fn single_sequence_canonical_kmers(filepath: String, k: usize) {}```  
//! Returns k-mer counts for individual sequences in a fasta file.  

use bio::{alphabets::dna::revcomp, io::fasta};
use dashmap::DashMap;
use fxhash::FxHasher;
use rayon::prelude::*;
use std::{
    env,
    error::Error,
    hash::BuildHasherDefault,
    io::{BufWriter, Write},
    str,
};

/// A simple struct for parsing command line k-size and filepath arguments.
pub struct Config {
    pub kmer_len: usize,
    pub filepath: String,
}

impl Config {
    pub fn new(mut args: env::Args) -> Result<Config, Box<dyn Error>> {
        args.next();

        let mut kmer_len: usize = 0;

        if let Some(arg) = args.next() {
            kmer_len = arg.as_str().parse()?;
        }

        let filepath = match args.next() {
            Some(arg) => arg,
            None => return Err("Issue with filepath input".into()),
        };
        Ok(Config { kmer_len, filepath })
    }
}

/// A custom `DashMap` w/ `FxHasher`.  
///  
/// ```use dashmap::DashMap;```  
/// ```use fxhash::FxHasher;```  
/// ```// skip```  
/// ```let dashfx_hash: DashFx = DashMap::with_hasher(BuildHasherDefault::<FxHasher>::default());```
/// Useful: [Using a Custom Hash Function in Rust](https://docs.rs/hashers/1.0.1/hashers/#using-a-custom-hash-function-in-rust).
pub type DashFx = DashMap<Box<[u8]>, u32, BuildHasherDefault<FxHasher>>;

///  - Reads sequences from fasta records in parallel using [`rayon`](https://docs.rs/rayon/1.5.1/rayon/),
/// using a customized [`dashmap`](https://docs.rs/dashmap/4.0.2/dashmap/struct.DashMap.html)
/// with [`FxHasher`](https://docs.rs/fxhash/0.2.1/fxhash/struct.FxHasher.html) to update in parallel a
/// hashmap of canonical k-mers (keys) and their frequency in the data (values).  
///  - Ignores substrings containing `N`.  
///  - Canonicalizes by lexicographically smaller of k-mer/reverse-complement.  
pub fn canonicalize_kmers(filepath: String, k: usize) -> Result<(), Box<dyn Error>> {
    let kmer_map: DashFx = DashMap::with_hasher(BuildHasherDefault::<FxHasher>::default());

    let _ = fasta::Reader::from_file(&filepath)?
        .records()
        .into_iter()
        .par_bridge()
        .for_each(|record| {
            let seq: &[u8] = record
                .as_ref()
                .unwrap_or_else(|_| {
                    panic!("{}", "problem getting fasta records from file".to_string())
                })
                .seq();

            for i in 0..(seq.len() + 1).saturating_sub(k) {
                match &seq[i..i + k] {
                    kmer if !kmer.contains(&78_u8) => match kmer {
                        kmer if kmer_map.contains_key(&Box::from(kmer)) => {
                            *kmer_map.get_mut(&Box::from(kmer)).unwrap() += 1
                        }

                        kmer if revcomp(kmer).as_slice() > kmer => {
                            *kmer_map.entry(Box::from(kmer)).or_insert(0) += 1
                        }

                        not_canonical if revcomp(not_canonical).as_slice() < not_canonical => {
                            *kmer_map
                                .entry(Box::from(revcomp(not_canonical)))
                                .or_insert(0) += 1
                        }

                        invalid if invalid.contains(&78_u8) => continue,

                        &_ => panic!("{}", "problem matching canonical kmer".to_string()),
                    },
                    &_ => continue,
                }
            }
        });

    let mut buf = BufWriter::new(std::io::stdout());

    kmer_map.into_iter().for_each(|(kmer, count)| {
        writeln!(buf, ">{}\n{}", count, str::from_utf8(&kmer).unwrap())
            .expect("Unable to write output");
    });
    buf.flush()?;
    Ok(())
}
