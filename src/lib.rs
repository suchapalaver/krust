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

pub fn run(config: Config) -> Result<(), Box<dyn Error>> {
    //  Get filepath and k-mer length
    let filepath: String = config.filepath;

    let k: usize = config.kmer_len;

    //  Create a hashmap of canonical k-mers and their frequency across all records
    let fasta_hash: DashMap<Box<[u8]>, u32> = kanonicalize(filepath, k);

    //  Print
    output(fasta_hash);

    Ok(())
}

pub fn kanonicalize(filepath: String, k: usize) -> DashMap<Box<[u8]>, u32> {
    let fasta_hash: DashMap<Box<[u8]>, u32> = DashMap::new();

    //  Read fasta records and count kmers
    fasta::Reader::from_file(&filepath).unwrap()
        .records()
        .into_iter()
        .par_bridge()
        .for_each(|record| {
            let seq: &[u8] = record.as_ref().unwrap().seq();

            for i in 0..(seq.len() + 1).saturating_sub(k) {
                //  Irradicate kmers containing 'N'
                if seq[i..i + k].contains(&b'N') {
                } else {
                    //  Canonicalize by lexicographically smaller of kmer/reverse-complement
                    let canon: Box<[u8]> =
                        Box::from(min(&seq[i..i + k], &revcomp(&seq[i..i + k])));
                    // update DashMap with canonical kmer count
                    *fasta_hash.entry(canon).or_insert(0) += 1;
                }
            }
        });
    fasta_hash
}

pub fn output(fasta_hash: DashMap<Box<[u8]>, u32>) {
    //  Create handle and BufWriter for writing
    let handle = std::io::stdout();

    let mut buf = BufWriter::new(handle);

    fasta_hash.into_iter().for_each(|(kmer, f)| {
        //  Write:
        //  >frequency across fasta file for both kmer and its reverse complement
        //  canonical k-mer
        writeln!(buf, ">{}\n{}", f, str::from_utf8(&kmer).unwrap())
            .expect("Unable to write data");
    });

    buf.flush().unwrap();
}
