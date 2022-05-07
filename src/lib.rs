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
//! - ```fn single_sequence_canonical_kmers(filepath: String, k: usize) {}```  
//! Returns k-mer counts for individual sequences in a fasta file.
//! - Testing!

use bio::io::fasta;
use dashmap::DashMap;
use fxhash::FxHasher;
use rayon::prelude::*;
use std::{
    env,
    error::Error,
    hash::BuildHasherDefault,
    io::{BufWriter, Write},
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
                Ok(kmer_len) if kmer_len > 0 && kmer_len < 33 => kmer_len,
                Ok(_) => return Err("k-mer length needs to be larger than zero and, for `krust` in its current working form, no more than 32".into()),
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

fn reverse(dna: &[u8]) -> Vec<u8> {
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

struct BP(u64);

impl BP {
    fn new(sub: &[u8]) -> BP {
	let bitpacked_kmer = sub.iter().fold(0, |mut k, byte| {
            k <<= 2;
	    let mask = match byte {
                b'A' => 0,
                b'C' => 1,
                b'G' => 2,
                b'T' => 3,
                _    => panic!("Won't happen!"),
            };
            k | mask
        });
	BP(bitpacked_kmer)
    }
}

struct PK((Vec<u8>, i32));

impl PK {
    fn new(pair: (u64, i32), k: usize) -> PK {
	let mut byte_string = Vec::new();
        for i in 0..k {
            let c = {
                let isolate = pair.0 << ((i * 2) + 64 - (k * 2));
                let base = isolate >> 62;
                match base {
                    0 => b'A',
                    1 => b'C',
                    2 => b'G',
                    3 => b'T',
                    _ => panic!("Won't ever happen!"),
                }
            };
            byte_string.push(c);
        }
        PK((byte_string, pair.1))
    }
}

/// A custom `DashMap` w/ `FxHasher`.  
///  
/// ```use dashmap::DashMap;```  
/// ```use fxhash::FxHasher;```  
/// ```// skip```  
/// ```let dashfx_hash: DashFx = DashMap::with_hasher(BuildHasherDefault::<FxHasher>::default());```
/// Useful: [Using a Custom Hash Function in Rust](https://docs.rs/hashers/1.0.1/hashers/#using-a-custom-hash-function-in-rust).
pub type DashFx = DashMap<u64, i32, BuildHasherDefault<FxHasher>>;

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

	    let mut i = 0;
	    while i <= seq.len() - k {
		let sub = &seq[i..i + k];
		if sub.contains(&b'N') {
		    i += k - 1;
		} else {
		    // get revcomp
		    let x = reverse(sub);
		    // bitpack 'canonical kmer'
		    let bitpacked_kmer = {
			if x.as_slice() < sub {
			    BP::new(&x)
			} else {
			    BP::new(sub)
			}};
		    *kmer_map.entry(bitpacked_kmer.0).or_insert(0) += 1;
		    i += 1;
		}
	    }
	});

    let mut buf = BufWriter::new(std::io::stdout());

    kmer_map.into_iter()
	.map(|pair| PK::new(pair, k).0)
        .for_each(|(kmer, count)| {
            writeln!(buf, ">{}\n{}", count, std::str::from_utf8(&kmer).unwrap())
        .expect("unable to write output");
    });
    buf.flush()?;
    
    Ok(())
}
