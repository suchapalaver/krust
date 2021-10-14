use bio::{alphabets::dna::revcomp, io::fasta};
use dashmap::DashMap;
use rayon::prelude::*;
use std::{
    cmp::min,
    env,
    error::Error,
    io::{BufWriter, Write},
    str,
};

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

///  Records canonical k-mers and their frequency across all records.
///  Uses DashMap to store canonical k-mers and their frequency in the data.
///  Reads sequences from fasta records in parallel using rayon (crate).
///  Ignores substrings containing 'N'.
///  Canonicalizes by lexicographically smaller of k-mer/reverse-complement
///  Prints to stdout.
///  Creates handle and BufWriter, writes on alternate lines,
///  ">{frequency across fasta file for both canonical k-mer and its reverse complement}",
///  "{canonical k-mer}"
pub fn run(config: Config) -> Result<(), Box<dyn Error>> {
    let k = config.kmer_len;
    assert!(
        k > 0,
        "The requested k-mer length caused a problem. k = {}",
        k
    );

    {
        let fasta_hash = canonicalize_kmers(config.filepath, k);
        let mut buf = BufWriter::new(std::io::stdout());

        fasta_hash.into_iter().for_each(|(kmer, count)| {
            writeln!(buf, ">{}\n{}", count, str::from_utf8(&kmer).unwrap())
                .expect("Unable to write output");
        });

        buf.flush().unwrap();
    };

    Ok(())
}

///  Reads sequences from fasta records in parallel using rayon (crate).
///  Ignores substrings containing 'N'.
///  Canonicalizes by lexicographically smaller of k-mer/reverse-complement
///  Returns a DashMap canonical k-mer keys and their respective counts in the data.
pub fn canonicalize_kmers(filepath: String, k: usize) -> DashMap<Box<[u8]>, u64> {
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

    canonical_hash
}

/*
Future:
/// Returns k-mer counts for individual sequences in a fasta file
pub fn single_sequence_canonical_kmers(filepath: String, k: usize) {}
*/
