use std::{
    env,
    error::Error,
    fs::{self, File},
    io::Write,
    path::Path,
    str
};
extern crate bio;
use bio::{alignment::sparse::hash_kmers,
	  alphabets::dna::revcomp,
	  io::fasta
};
extern crate rayon;
use rayon::{iter::ParallelBridge,
	    prelude::ParallelIterator
};
use std::time::{Duration,
		Instant
};

pub struct Config {
    pub kmer_len: String,
    pub filepath: String,
}

impl Config {
    
    pub fn new(mut args: env::Args) -> Result<Config, &'static str> {

	args.next();

        let kmer_len = match args.next() {
            Some(arg) => arg,
            None => return Err("Didn't get a k-mer length"),
        };
        let filepath = match args.next() {
            Some(arg) => arg,
            None => return Err("Didn't get a file name"),
        };
        Ok(Config { kmer_len, filepath })
    }
}

pub fn run(config: Config) -> Result<(), Box<dyn Error>> {
    let start = Instant::now();
    
    fs::create_dir("output")?;

    let filepath: String = config.filepath;

    let kmer_len = config.kmer_len.parse::<usize>().unwrap();

    let reader = fasta::Reader::from_file(&filepath).unwrap();

    reader.records().par_bridge().for_each(|result| {
	
        let result_data = result.as_ref().unwrap();

        let pathname = format!("output/{}.tsv", result_data.id());

        let path = Path::new(&pathname);

        let display = path.display();

        let mut file = match File::create(&path) {
	    
            Err(why) => panic!("couldn't create {}: {}", display, why),
            Ok(file) => file,
        };
	
        for (kmer, kmer_positions) in hash_kmers(result_data.seq(), kmer_len) {

	    let rvc = revcomp(kmer);

            match str::from_utf8(kmer) {
		
                Err(e) => println!("Problem: {}", e),
                Ok(kmer_s) => match str::from_utf8(&rvc) {

		    Err(why) => panic!("couldn't write to {}: {}", display, why),
		    Ok(rvc) => {

			let data = format!("{}\t{}\t{}\n", kmer_s, rvc, kmer_positions.len());
                        write!(file, "{}", data).expect("Unable to write file");
                    }
                },
            }
        }
    });
    let duration = start.elapsed();
    println!("Time elapsed in expensive_function() is: {:?}", duration);
    Ok(println!("{}", filepath))
}
