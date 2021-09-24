use bio::{alphabets::dna::revcomp, io::fasta};
use dashmap::DashMap;
use rayon::prelude::*;
use std::{env, error::Error, fs::File, io::Write, str, time::Instant};
use std::io::BufWriter;

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
    //  Data filepath
    let filepath: String = config.filepath;
    //  K-mer length k
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
    fasta_records.par_iter().for_each(|result| {
        let result_data: &fasta::Record = result.as_ref().unwrap();
        {
            let seq = result_data.seq();

            for i in 0..(seq.len() + 1).saturating_sub(k) {
                fasta_hash
                    .entry(&seq[i..i + k])
                    .or_insert_with(Vec::new)
                    .push(i as u32);
            }
        }
    });
    //  End of benchmark timer
    let uniq_duration = hash_start.elapsed();
    eprintln!(
        "Time elapsed merging sequence kmer data into hashmap: {:?}\n",
        uniq_duration
    );

    //  PRINTING OUTPUT

    //  Benchmark timer
    let print_start = Instant::now();

    //  Create handle for writing
    //let handle = File::create("output.tsv").unwrap();

    let handle = &std::io::stdout();

    let mut buf = BufWriter::new(handle);
    

    fasta_hash.into_iter().for_each(|(k, f)| {
        //  Convert k-mer bytes to str
        let kmer = str::from_utf8(k).unwrap();

        //  Don't write k-mers containing 'N'
        if k.contains(&b'N') {
        } else {
            //  Use bio (crate) revcomp to get k-mer reverse complement
            let rvc: Vec<u8> = revcomp(k as &[u8]);

            //  Convert revcomp from bytes to str
            let rvc: &str = str::from_utf8(&rvc).unwrap();

            //  Write (separated by tabs):
            //        k-mer
            //        reverse complement
            //        frequency across fasta file
	    writeln!(buf, "{}\t{}\t{}", kmer, rvc, f.len()).expect("Unable to write data");
        }
    });
    buf.flush().unwrap();
    //  END OF WRITING OUTPUT
    let duration = print_start.elapsed();

    eprintln!("Time elapsed printing output: {:?}\n", duration);

    Ok(())
}
