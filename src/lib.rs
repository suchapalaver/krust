use bio::{alignment::sparse::hash_kmers, alphabets::dna::revcomp, io::fasta};
use rayon::prelude::*;
use std::{
    env,
    error::Error,
    fs::{self, File},
    io::Write,
    path::Path,
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
        let filepath = match args.next() {
            Some(arg) => arg,
            None => return Err("Didn't get a file name"),
        };
        Ok(Config { kmer_len, filepath })
    }
}

pub fn run(config: Config) -> Result<(), Box<dyn Error>> {
    fs::create_dir("output")?;

    let filepath: String = config.filepath;

    let kmer_len: usize = config.kmer_len;

    let reader: fasta::Reader<std::io::BufReader<File>> =
        fasta::Reader::from_file(&filepath).unwrap();

    let fasta_records: Vec<Result<fasta::Record, std::io::Error>> = reader.records().collect();

    fasta_records.par_iter().for_each(|result| {
        let result_data: &fasta::Record = result.as_ref().unwrap();

        let pathname: String = format!("output/{}.tsv", result_data.id());

        let path: &Path = Path::new(&pathname);

        let mut file: File = File::create(&path).expect("Couldn't create a file");

        for (kmer, kmer_positions) in hash_kmers(result_data.seq(), kmer_len) {
            let kmer_s: String = kmer.iter().map(|c| *c as char).collect::<String>();

	    let rvc: String = revcomp(kmer).iter().map(|c| *c as char).collect::<String>();

            write!(file, "{}\t{}\t{}\n", kmer_s, rvc, kmer_positions.len())
                .expect("Unable to write file");
        }
    });
    Ok(println!("{}", filepath))
}
