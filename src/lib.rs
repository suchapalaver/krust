use bio::{alphabets::dna::revcomp, io::fasta};
use dashmap::DashMap;
use rayon::prelude::*;
use std::{
    cmp::min,
    env,
    error::Error,
    fs::File,
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
    //  Get parameters
    //  Fasta filepath
    let filepath: String = config.filepath;
    //  K-mer length k
    let k: usize = config.kmer_len;

    //  Create fasta file reader
    let reader: fasta::Reader<std::io::BufReader<File>> =
        fasta::Reader::from_file(&filepath).unwrap();

    //  Read fasta records into a vector
    let fasta_records: Vec<Result<fasta::Record, std::io::Error>> = reader.records().into_iter().par_bridge().collect();

    eprintln!("Number of records in fasta file: {}\n", fasta_records.len());

    //  Create a Dashmap, a hashmap mutably accessible from different parallel processes
    let fasta_hash: DashMap<Vec<u8>, Vec<u32>> = DashMap::new();

    //  Iterate through fasta records in parallel
    fasta_records.par_iter().for_each(|result| {
        let result_data: &fasta::Record = result.as_ref().unwrap();

        let seq: &[u8] = result_data.seq();

        for i in 0..(seq.len() + 1).saturating_sub(k) {
            let kmer: Vec<u8> = min(seq[i..i + k].to_vec(), revcomp(&seq[i..i + k]));

            fasta_hash
                .entry(kmer)
                .or_insert_with(Vec::new)
                .push(i as u32);
        }
    });

    //  PRINTING OUTPUT

    //  Create handle and BufWriter for writing
    let handle = &std::io::stdout();

    let mut buf = BufWriter::new(handle);

    for (k, f) in fasta_hash.into_iter() {
        //  Convert k-mer bytes to str
        let kmer = str::from_utf8(&k).unwrap();

        //  Don't write k-mers containing 'N'
        if kmer.contains('N') {
        } else {
            //  Write (separated by tabs):
            //  k-mer (lexicographically smaller of k-mer, reverse complement pair)
            //  frequency across fasta file for both kmer and its reverse complement
            writeln!(buf, ">{}\n{}", f.len(), kmer).expect("Unable to write data");
        }
    }
    buf.flush().unwrap();

    Ok(())
}
