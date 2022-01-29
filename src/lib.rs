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

use bio::io::fasta;
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

/// Parsing command line k-size and filepath arguments.
pub struct Config {
    pub kmer_len: usize,
    pub filepath: String,
}

impl Config {
    pub fn new(mut args: env::Args) -> Result<Config, Box<dyn Error>> {
        let kmer_len: usize = match args.nth(1) {
            Some(arg) => match arg.parse() {
                Ok(kmer_len) if kmer_len > 0 => kmer_len,
                Ok(_) => return Err("k-mer length needs to be larger than zero".into()),
                Err(_) => return Err(format!("issue with k-mer length argument: {}", arg).into()),
            },
            None => return Err("k-mer length input required".into()),
        };

        let filepath = match args.next() {
            Some(arg) => arg,
            None => return Err("filepath argument needed".into()),
        };

        Ok(Config { kmer_len, filepath })
    }
}

// to do:
/*
pub fn run(filepath: String, k: usize) -> Result<(), Box<dyn Error>> {
    let kmer_map: DashFx = DashMap::with_hasher(BuildHasherDefault::<FxHasher>::default());

    let _ = fasta::Reader::from_file(&filepath)?
        .records()
        .into_iter()
        .par_bridge()
        .for_each(|r| {
            let record = r.expect("error reading fasta record");

            let seq: &[u8] = record.seq();

        canonicalize_kmers(seq, k)?;
    });
    Ok(())
}

fn canonicalize_kmers(seq: &[u8], k: usize) -> Result<(), Box<dyn Error>> {

}
 */
const VALID_BYTES: &[u8] = &[65_u8, 67_u8, 71_u8, 84_u8];

fn dna_strand(dna: &[u8]) -> Vec<u8> {
    let revcomp = dna
        .iter()
        .rev()
        .map(|c| match *c {
            67_u8 => 71_u8,
            71_u8 => 67_u8,
            84_u8 => 65_u8,
            _ => 84_u8, //65_u8
        })
        .collect();
    revcomp
}

fn is(bytes: &[u8]) -> bool {
    bytes
        .iter()
	.all(|x| VALID_BYTES.contains(x))
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
        .for_each(|r| {
            let record = r.expect("error reading fasta record");

            let seq: &[u8] = record.seq();

            for i in 0..(seq.len() + 1).saturating_sub(k) {
                if kmer_map.contains_key(&Box::from(&seq[i..i + k])) {
                    *kmer_map.get_mut(&Box::from(&seq[i..i + k])).unwrap() += 1;
                    continue;
                }

                match &seq[i..i + k] {
                    valid if is(valid) => match dna_strand(valid) {
                        x if x.as_slice() < &seq[i..i + k] => {
                            *kmer_map.entry(Box::from(x)).or_insert(0) += 1
                        }
                        _ => {
   
                            *kmer_map.entry(Box::from(&seq[i..i + k])).or_insert(0) += 1
                        }
                    },
                    _ => continue,
                }
            }
        });

    let mut buf = BufWriter::new(std::io::stdout());

    kmer_map.into_iter().for_each(|(kmer, count)| {
        writeln!(
            buf,
            ">{}\n{}",
            count,
            str::from_utf8(&kmer).expect("couldn't convert k-mer to readable format")
        )
        .expect("unable to write output");
    });
    buf.flush()?;
    Ok(())
}
