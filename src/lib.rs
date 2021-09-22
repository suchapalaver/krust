use bio::{alignment::sparse::hash_kmers, alphabets::dna::revcomp, io::fasta};
use dashmap::DashMap;
use rayon::prelude::*;
use std::{env, error::Error, fs::File, io::Write, str, time::Instant};
use fxhash::{FxHashMap, FxHashSet};

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

pub fn hash_fasta_rec(
    result: &Result<fasta::Record, std::io::Error>,
    k: usize,
) -> FxHashMap<&[u8], usize> {
    let result_data: &fasta::Record = result.as_ref().unwrap();

    let mut new_hashmap = FxHashMap::default();

    for (kmer, kmer_pos) in hash_kmers(result_data.seq(), k) {
        new_hashmap.insert(kmer, kmer_pos.len());
    }
    new_hashmap
}

pub fn run(config: Config) -> Result<(), Box<dyn Error>> {
    let start = Instant::now();

    let filepath: String = config.filepath;

    let k: usize = config.kmer_len;

    let reader: fasta::Reader<std::io::BufReader<File>> =
        fasta::Reader::from_file(&filepath).unwrap();

    let fasta_records: Vec<Result<fasta::Record, std::io::Error>> = reader.records().collect();

    let mut hash_vec: Vec<FxHashMap<&[u8], usize>> = fasta_records
        .par_iter()
        .map(|result| hash_fasta_rec(result, k))
        .collect();
    
    let hash_duration = start.elapsed();

    // merging hashmaps

    let mut hash_len_vec = FxHashSet::default(); // create set of number of kmers 
    
    for h in &hash_vec {
	hash_len_vec.insert(h.len());
    }
    let longest_len = hash_len_vec.iter().max().unwrap();
    
    let i = &hash_vec.iter().position(|h| h.len() == *longest_len).unwrap();

    let this = hash_vec.remove(*i);
    
    let final_hash: DashMap<&[u8], usize> = DashMap::default();

    this.par_iter().for_each(|(k, v)| { final_hash.insert(k, *v); });

    hash_vec.par_iter().for_each(|h| {
        for (kmer, freq) in h {
            if final_hash.contains_key(kmer) {
		*final_hash.get_mut(kmer).unwrap() += freq;
            } else {
                final_hash.insert(kmer, *freq);
            }
        }
    });

    let uniq_duration = start.elapsed();

    let stdout_ref = &std::io::stdout();

    final_hash.into_iter().par_bridge().for_each(|(k, f)| {
        let kmer = str::from_utf8(k).unwrap();

        if kmer.contains("N") {
        } else {
            let rvc = revcomp(k);

            let rvc = str::from_utf8(&rvc).unwrap();

            let mut lck = stdout_ref.lock();

            writeln!(&mut lck, "{}\t{}\t{}", kmer, rvc, f).expect("Couldn't write output");
        }
    });
    let duration = start.elapsed();

    eprintln!(
        "Time elapsed creating hashmaps of all kmers in all sequences: {:?}\n",
        hash_duration
    );
    eprintln!("Time elapsed merging hashmaps: {:?}\n", uniq_duration);

    eprintln!("Time elapsed in runtime: {:?}\n", duration);

    Ok(())
}
