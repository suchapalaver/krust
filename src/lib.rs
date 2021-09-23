use bio::{alphabets::dna::revcomp, io::fasta};
use dashmap::DashMap;
use rayon::{iter::ParallelBridge, prelude::*};
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
    //  Get parameters
    let filepath: String = config.filepath;

    let k: usize = config.kmer_len;

    //  Create fasta file reader
    let reader: fasta::Reader<std::io::BufReader<File>> =
        fasta::Reader::from_file(&filepath).unwrap();

    //  Read fasta records into a vector
    let fasta_records: Vec<Result<fasta::Record, std::io::Error>> = reader.records().collect();

    eprintln!("Number of records in fasta file: {}\n", fasta_records.len());

    //  Benchmark timer
    let hash_start = Instant::now();
    
    //  Create a Dashmap, a hashmap mutably accessible from different parallel processes
    let fasta_hash: DashMap<&[u8], Vec<u32>> = DashMap::new();
    
    //  Iterate through fasta records in parallel
    fasta_records
        .par_iter() 
        .for_each(|result| {
            let result_data: &fasta::Record = result.as_ref().unwrap();
            {
                let seq = result_data.seq();

		let slc = seq;
                
                for i in 0..(slc.len() + 1).saturating_sub(k) {
                    fasta_hash.entry(&slc[i..i + k])
                        .or_insert_with(Vec::new)
                        .push(i as u32);
                }
            }
        });
    //  End of benchmark timer
    let uniq_duration = hash_start.elapsed();
    eprintln!("Time elapsed merging sequence kmer data into hashmap: {:?}\n", uniq_duration);

    //  PRINTING OUTPUT -- by far the slowest task

    //  Benchmark timer
    let print_start = Instant::now();

    //  Create handle for writing to standard output
    let stdout_ref = &std::io::stdout();
    
    // Iterate in parallel through fasta_hash with rayon (crate) parallel bridge 
    fasta_hash.into_iter().par_bridge().for_each(|(k, f)| {
	//  Convert kmer bytes to str
        let kmer = str::from_utf8(k).unwrap();

	//  Don't write kamers containing 'N'
        if kmer.contains('N') {
        } else {
	    //  Use bio (crate) revcomp to get kmer reverse complement
            let rvc: Vec<u8> = revcomp(k as &[u8]);

	    //  Convert revcomp from bytes to str
            let rvc = str::from_utf8(&rvc).unwrap();

	    //  Create mutable lock to write to standard output from paralell process
            let mut lck = stdout_ref.lock();

	    //  Write (separated by tabs):
	    //        kmer
	    //        reverse complement
	    //        frequency across fasta file 
            writeln!(&mut lck, "{}\t{}\t{}", kmer, rvc, f.len()).expect("Couldn't write output");
        }
    });
    //  END OF WRITING OUTPUT
    let duration = print_start.elapsed();

    eprintln!("Time elapsed printing output: {:?}\n", duration);

    Ok(())
}
