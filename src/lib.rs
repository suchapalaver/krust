use bio::{alphabets::dna::revcomp, io::fasta};
use dashmap::DashMap;
use rayon::prelude::*;
use std::{cmp::{max, min}, env, error::Error, fs::File, io::{BufWriter, Write}, rc::Rc, str, time::Instant};


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
/*
struct RevComps {
    revcomps: DashMap<Vec<u8>, Vec<u8>>
}

impl RevComps {
    pub fn new() -> Self {
        RevComps { revcomps: DashMap::new() }
    }
    
    pub fn insert(&mut self, kmer: Vec<u8>, revcomp: Vec<u8>) {
        self.revcomps.insert(kmer, revcomp);
    }
    /*
    pub fn print_all(&self) {
        self.revcomps.iter().for_each(|item| println!("{:#?}", item));
    }
     */
}
*/
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
    let fasta_hash: DashMap<Vec<u8>, Vec<u32>> = DashMap::new();

    let revcomp_hash: DashMap<Vec<u8>, Vec<u8>> = DashMap::new();

    //  Iterate through fasta records in parallel
    fasta_records.par_iter().for_each(|result| {
        
        {
	    let result_data: &fasta::Record = result.as_ref().unwrap();
	    
            let seq = result_data.seq();
	    
	    for i in 0..(seq.len() + 1).saturating_sub(k) {

		let kmer = min(seq[i..i + k].to_vec(), revcomp(&seq[i..i + k]));

		let rvc: Vec<u8> = max(seq[i..i + k].to_vec(), revcomp(&seq[i..i + k]));
		
		if revcomp_hash.contains_key(&kmer) {
		} else {
		    revcomp_hash.insert(kmer, rvc);
		}
	
		fasta_hash
		    .entry(seq[i..i + k].to_vec())
		    .or_insert_with(Vec::new)
		    .push(i as u32);
	
	    }
		//let kmer = min(&seq[i..i + k], rc_rv);  //  idea: identify kmer as opp to revcomp by lex.comp., then only update revcomp 
	}
    });
    //  End of benchmark timer
    let uniq_duration = hash_start.elapsed();
    eprintln!(
        "Time elapsed merging sequence kmer data into hashmap: {:?}\n",
        uniq_duration
    );
/*
    let final_hash = fasta_hash.into_iter().for_each(|(k, f)| {
	let rvc = revcomp(k);
	let rc = &rvc;
	let kmer = min(k, &rc);
	let revcomp = max(k, &rc);
	
	let rev_vals = fasta_hash.remove(&revcomp).unwrap(); // add this to vals for kmer
	f.extend(rev_vals.1);
	
    });
    */
    //  PRINTING OUTPUT

    //  Benchmark timer
    let print_start = Instant::now();

    //  Create handle and BufWriter for writing
    let handle = &std::io::stdout();

    let mut buf = BufWriter::new(handle);

    fasta_hash.into_iter().for_each(|(k, f)| {
        //  Convert k-mer bytes to str
        let kmer = str::from_utf8(&k).unwrap();
	
        //  Don't write k-mers containing 'N'
        if kmer.contains('N') {
        } else {
	    //  Use bio (crate) revcomp to get k-mer reverse complement
	    //let rvc: Vec<u8> = revcomp(k as &[u8]);
	    
	    //  Convert revcomp from bytes to str
	    //let rvc: &str = str::from_utf8(&rvc).unwrap();
	    
	    //  Write (separated by tabs):
	    //        k-mer
	    //        reverse complement
	    //        frequency across fasta file
	    //let rv_vals = fasta_hash.get(revcomp).unwrap().len();	
	    writeln!(buf, "{}\t{}", kmer, f.len()).expect("Unable to write data");
	}	
    });
    buf.flush().unwrap();


    //  Create handle and BufWriter for writing
    let handle = &std::io::stdout();

    let mut buf = BufWriter::new(handle);

    revcomp_hash.into_iter().for_each(|(k, r)| {
        //  Convert k-mer bytes to str
        let kmer = str::from_utf8(&k).unwrap();
	let rvc = str::from_utf8(&r).unwrap();
	
        //  Don't write k-mers containing 'N'
        if kmer.contains('N') {
        } else {
	    //  Use bio (crate) revcomp to get k-mer reverse complement
	    //let rvc: Vec<u8> = revcomp(k as &[u8]);
	    
	    //  Convert revcomp from bytes to str
	    //let rvc: &str = str::from_utf8(&rvc).unwrap();
	    
	    //  Write (separated by tabs):
	    //        k-mer
	    //        reverse complement
	    //        frequency across fasta file
	    //let rv_vals = fasta_hash.get(revcomp).unwrap().len();	
	    writeln!(buf, "{}\t{}", kmer, rvc).expect("Unable to write data");
	}	
    });
    buf.flush().unwrap();

    
    //  END OF WRITING OUTPUT

    let duration = print_start.elapsed();

    eprintln!("Time elapsed printing output: {:?}\n", duration);

    Ok(())
}
