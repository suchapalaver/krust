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
pub fn run(filepath: String, k: usize) -> Result<(), Box<dyn Error>> {
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
		let bytestring = Kmer::try_from(sub);
		match bytestring.is_ok() {
                    true => {
                        let Kmer(bytestring) = bytestring.unwrap();
                        let RevCompKmer(revcompkmer) = RevCompKmer::from(&bytestring);
                        let CanonicalKmer(canonical_kmer) =
                            CanonicalKmer::from((revcompkmer, bytestring));
                        let BitpackedKmer(kmer) = BitpackedKmer::from(canonical_kmer);
                        *kmer_map.entry(kmer).or_insert(0) += 1;
                        i += 1;
                    }
                    false => {
			let mut position: usize = 0;
			sub.iter().enumerate().for_each(|(pos, byte)| {
			    if byte != &(b'A' | b'C' | b'G' | b'T') {
				position = pos;
			    }
			});
			i += position + 1;
		    },
		}
	    }
        });

    let mut buf = BufWriter::new(std::io::stdout());
    kmer_map
        .into_iter()
        .map(|(kmer, freq)| (UnpackedKmer::from((kmer, k)), freq))
        .for_each(|(UnpackedKmer(kmer), count)| {
            writeln!(
                buf,
                ">{}\n{}",
                count,
                std::str::from_utf8(kmer.as_slice()).unwrap()
            )
            .expect("Unable to write output.");
        });
    buf.flush()?;

    Ok(())
}
#[derive(Debug)]
pub struct Kmer(Vec<u8>);

impl TryFrom<&[u8]> for Kmer {
    type Error = ();
    fn try_from(sub: &[u8]) -> Result<Kmer, ()> {
        match !sub.contains(&b'N') {
            true => Ok(Kmer(sub.to_vec())),
            false => Err(()),
        }
    }
}

/// Compressing k-mers of length `0 < k < 33`, bitpacking them into unsigned integers.
pub struct BitpackedKmer(u64);

impl From<Vec<u8>> for BitpackedKmer {
    fn from(sub: Vec<u8>) -> Self {
        let bitpacked_kmer: u64 = {
            let mut k: u64 = 0;
            for byte in sub.iter() {
                k <<= 2;
                let mask = match *byte {
                    b'A' => 0,
                    b'C' => 1,
                    b'G' => 2,
                    _ => 3, // b'T'
                };
                k |= mask;
            }
            k
        };
        BitpackedKmer(bitpacked_kmer)
    }
}

/// Converting a DNA string slice into its [reverse compliment](https://en.wikipedia.org/wiki/Complementarity_(molecular_biology)#DNA_and_RNA_base_pair_complementarity).
pub struct RevCompKmer(Vec<u8>);

impl From<&Vec<u8>> for RevCompKmer {
    fn from(sub: &Vec<u8>) -> Self {
        RevCompKmer(
            sub.iter()
                .rev()
                .map(|byte| RevCompKmer::complement(*byte))
                .collect::<Vec<u8>>(),
        )
    }
}

impl RevCompKmer {
    fn complement(byte: u8) -> u8 {
        match byte {
            67_u8 => 71_u8,
            71_u8 => 67_u8,
            84_u8 => 65_u8,
            _ => 84_u8, // 65_u8
        }
    }
}

pub struct CanonicalKmer(Vec<u8>);

impl From<(Vec<u8>, Vec<u8>)> for CanonicalKmer {
    fn from(comp: (Vec<u8>, Vec<u8>)) -> Self {
        let canonical_kmer = match comp.0 < comp.1.to_vec() {
            true => comp.0,
            false => comp.1.to_vec(),
        };
        CanonicalKmer(canonical_kmer)
    }
}

/// Unpacking compressed, bitpacked k-mer data.
pub struct UnpackedKmer(Vec<u8>);

impl From<(u64, usize)> for UnpackedKmer {
    fn from(kmer_data: (u64, usize)) -> Self {
        let mut byte_string = Vec::new();
        let (kmer, k) = (kmer_data.0, kmer_data.1);
        for i in 0..k {
            let isolate = kmer << ((i * 2) + 64 - (k * 2));
            let base = isolate >> 62;
            let byte = UnpackedKmerByte::from(base);
            byte_string.push(byte.0);
        }
        UnpackedKmer(byte_string)
    }
}

/// Unpacking compressed, bitpacked k-mer data.
pub struct UnpackedKmerByte(u8);

impl From<u64> for UnpackedKmerByte {
    fn from(base: u64) -> Self {
        match base {
            0 => UnpackedKmerByte(b'A'),
            1 => UnpackedKmerByte(b'C'),
            2 => UnpackedKmerByte(b'G'),
            _ => UnpackedKmerByte(b'T'), // 3
        }
    }
}

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
