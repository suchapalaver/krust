use bio::{alignment::sparse::hash_kmers, alphabets::dna::revcomp, io::fasta};
use dashmap::DashMap;
use rayon::iter::ParallelBridge;
use rayon::prelude::*;
use std::{env, error::Error, fs::File, io::Write, str, time::Instant};

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
    let start = Instant::now();

    let filepath: String = config.filepath;

    let k: usize = config.kmer_len;

    let reader: fasta::Reader<std::io::BufReader<File>> =
        fasta::Reader::from_file(&filepath).unwrap();

    let fasta_records: Vec<Result<fasta::Record, std::io::Error>> = reader.records().collect();

    eprintln!("Number of records in fasta file: {}\n", fasta_records.len());

    let final_hash: DashMap<&[u8], usize> = DashMap::new();

    fasta_records
        .par_iter() //  "Where you use par_iter(), instead of using map try using fold then reduce. With rayon, fold will let you merge the data into HashMaps in parallel. Then reduce will take those maps and let you merge them into a single one. If you want to avoid the Vec altogether, you should be able to call par_bridge directly on the records() result instead of calling collect (then call fold and reduce). par_bridge creates a parallel iterator from a regular iterator."
        .for_each(|result| {
            let result_data: &fasta::Record = result.as_ref().unwrap();

            for (kmer, kmer_pos) in hash_kmers(result_data.seq(), k) {
                // the Vec<u32> that hash_kmers (rust-bio) is a list of indices of kmer's positions in seq.
                final_hash
                    .entry(kmer)
                    .and_modify(|v| *v += kmer_pos.len())
                    .or_insert(kmer_pos.len());
            }
        });

    let uniq_duration = start.elapsed();

    eprintln!("Time elapsed merging hashmaps: {:?}\n", uniq_duration);

    // END OF MERGING

    let stdout_ref = &std::io::stdout();

    final_hash.into_iter().par_bridge().for_each(|(k, f)| {
        let kmer = str::from_utf8(k).unwrap();

        if kmer.contains('N') {
        } else {
            let rvc = revcomp(k as &[u8]);

            let rvc = str::from_utf8(&rvc).unwrap();

            let mut lck = stdout_ref.lock();

            writeln!(&mut lck, "{}\t{}\t{}", kmer, rvc, f).expect("Couldn't write output");
        }
    });

    let duration = start.elapsed();

    eprintln!("Time elapsed in runtime: {:?}\n", duration);

    Ok(())
}
