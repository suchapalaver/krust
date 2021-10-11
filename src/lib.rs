use bio::{alphabets::dna::revcomp, io::fasta};
use dashmap::{DashMap, DashSet};
use rayon::prelude::*;
use std::{
    cmp::{max, min},
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

pub fn run(config: Config) -> Result<(), Box<dyn Error>> {
    //  Get filepath and k-mer length
    //  Create a hashmap of canonical k-mers and their frequency across all records
    //  Print output

    let filepath: String = config.filepath;

    let k: usize = config.kmer_len;

    let fasta_hash: DashMap<Box<[u8]>, u32> = canonicalize(filepath, k);

    output(fasta_hash);

    Ok(())
}

pub fn canonicalize(filepath: String, k: usize) -> DashMap<Box<[u8]>, u32> {
    //  Use DashMap for storing canonical k-mers and their frequency in the data.
    //  Use DashSet to log reverse compliment k-mers.
    //  Read sequences from fasta records in parallel using rayon (crate).
    //  Ignore substrings containing 'N'.
    //  Canonicalize by lexicographically smaller of k-mer/reverse-complement
    //  Skip performing repeat lexicographical size / reverse complement conversions.

    let canonical_hash: DashMap<Box<[u8]>, u32> = DashMap::new();

    let lookup: DashSet<Box<[u8]>> = DashSet::new();

    fasta::Reader::from_file(&filepath)
        .unwrap()
        .records()
        .into_iter()
        .par_bridge()
        .for_each(|record| {
            let seq: &[u8] = record.as_ref().unwrap().seq();

            for i in 0..(seq.len() + 1).saturating_sub(k) {
                let substring: &[u8] = &seq[i..i + k];

                if substring.contains(&b'N') {
                } else {
                    let boxed_sub: Box<[u8]> = Box::from(substring);

                    if lookup.contains(&boxed_sub) {
                        //let canon: Box<[u8]>= Box::from(&**lookup.get(&substring).unwrap().value());
                        *canonical_hash.entry(Box::from(revcomp(substring))).or_insert(0) += 1;
                    } else if let Some(mut x) = canonical_hash.get_mut(&boxed_sub) {
			*x += 1;
		    } else {
                        *canonical_hash.entry(Box::from(min(substring, &revcomp(substring)))).or_insert(0) += 1;

                        lookup.insert(Box::from(max(substring, &revcomp(substring))));
                    }
                }
            }
        });
    canonical_hash
}

pub fn output(fasta_hash: DashMap<Box<[u8]>, u32>) {
    //  Create handle and BufWriter and write on alternate lines:
    //  ">{frequency across fasta file for both canonical k-mer and its reverse complement}"
    //  "{canonical k-mer}"
    let handle = std::io::stdout();

    let mut buf = BufWriter::new(handle);

    fasta_hash.into_iter().for_each(|(kmer, f)| {
        writeln!(buf, ">{}\n{}", f, str::from_utf8(&kmer).unwrap()).expect("Unable to write data");
    });

    buf.flush().unwrap();
}
