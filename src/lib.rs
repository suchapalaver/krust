use bio::{alignment::sparse::hash_kmers, alphabets::dna::revcomp, io::fasta};
use itertools::Itertools;
use rayon::prelude::*;
use rayon::iter::ParallelBridge;
use rayon::prelude::ParallelIterator;
use reduce::Reduce;
//use rustc_hash::FxHashMap; //  Check Stack Overflow for suggestions
use std::{collections::{HashMap, HashSet}, env, error::Error, fs::File, io::Write, str, time::Instant};

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
pub fn div_n_merge(hash_vec_a: Vec<HashMap<&[u8], usize>>, hash_vec_b: Vec<HashMap<&[u8], usize>>) -> Vec<HashMap<&[u8], usize>> {

    let mut a = Vec::new();
    
    let mut b = Vec::new();
	
    while hash_vec.len() > 1 {
	if a.is_empty() | b.is_empty() { // is | correct here?
	    
	    let mid = hash_vec.len()/2;
	    
	    for i in 0..=mid {
		
                let c = hash_vec[i].clone();
		
		a.push(c);
	    }
	    for i in mid..(hash_vec.len()+1) {
		let d = hash_vec[i].clone();
		
		b.push(d);
	    }
	} else {  //  if we've already done work on a and b ...
	    if a.len() > 1 { // do a first
		let mid = a.len()/2;
	    
		for i in 0..=mid {
		    let c = a[i].clone();
		
		    a.push(c);
	    }
	    for i in mid..(hash_vec.len()+1) {
		let d = hash_vec[i].clone();
		
		b.push(d);
	    }
*/
pub fn hash_fasta_rec(
    result: &Result<fasta::Record, std::io::Error>,
    k: usize,
) -> HashMap<&[u8], usize> {
    let result_data: &fasta::Record = result.as_ref().unwrap();

    let mut new_hashmap = HashMap::new();

    for (kmer, kmer_pos) in hash_kmers(result_data.seq(), k) { // rust-bio's hash_kmers function, returns iterator of tuples (&[u8], Vec<u32>), the Vec being a list of indices of positions of kmer. 
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
    
    /*
    let hash_vec: Vec<HashMap<&[u8], usize>> = fasta_records
        .par_iter() //  Where you use par_iter(), instead of using map try using fold then reduce. With rayon, fold will let you merge the data into HashMaps in parallel. Then reduce will take those maps and let you merge them into a single one. If you want to avoid the Vec altogether, you should be able to call par_bridge directly on the records() result instead of calling collect (then call fold and reduce). par_bridge creates a parallel iterator from a regular iterator.
        .map(|result| hash_fasta_rec(result, k))
        .collect();
*/
    let final_hash: HashMap<&[u8], usize> = fasta_records.par_iter()
	.map(|result| hash_fasta_rec(result, k)
	)
	.fold(||HashMap::new(), |mut a: HashMap<&[u8], usize>, b| {
	    a.extend(b.into_iter()); a}
	)
	.reduce(||HashMap::new(),|mut a, b| {
	    for (k, v) in b {
		if a.contains_key(k) {
		    a.insert(k, v + a[k]);
		} else {
		    a.insert(k, v);
		}
	    }
	    a}
	);
    
    let hash_duration = start.elapsed();
    let uniq_duration = start.elapsed();


    //eprintln!("number of hashmaps in vec: {}", hash_vec.len());
    
    // MERGING HASHMAPS
    /*
    // THIS WORKS
    hash_vec.into_iter().for_each(|h| {
        for (kmer, freq) in h {
            if final_hash.contains_key(kmer) {
                final_hash.insert(kmer, freq + final_hash[kmer]);
            } else {
                final_hash.insert(kmer, freq);
            }
        }
    });
     */
    /*
    //  THIS WORKS
    let final_hash: HashMap<&[u8], usize> = Reduce::reduce(hash_vec.into_iter(), |mut ha, hb| {
	    for (kmer, freq) in hb {
		if ha.contains_key(kmer) {
                    ha.insert(kmer, freq + ha[kmer]);
		} else {
                    ha.insert(kmer, freq);
		}
            }
	ha
    }).unwrap();
     */
    /*
    let final_hash: HashMap<&&[u8], usize> = hash_vec
	.par_iter()
	.fold(||HashMap::new(), |mut a, b| {
	    for (k, v) in b {
		if a.contains_key(k) {
		    eprintln!("cnt");
		    a.insert(k, v + a[k]);
		} else {
		    eprintln!("else");
		    a.insert(k, *v);
		}
	    }
	    //eprintln!("{:?}", &a);
	    //eprintln!("{:?}", &b);
	    a
	}).reduce(||HashMap::new(),
		  |mut a, b| {
		      for (k, v) in b {
			  if a.contains_key(k) {
			      a.insert(k, v + a[k]);
			  } else {
			      a.insert(k, v);
			  }
		      }
		      a
		  });
    */
    // END OF MERGING
    

    
/*
    let stdout_ref = &std::io::stdout();

    final_hash.par_iter().for_each(|(k, f)| {
        let kmer = str::from_utf8(k).unwrap();

        if kmer.contains("N") {
        } else {
            let rvc = revcomp(k as &[u8]);

            let rvc = str::from_utf8(&rvc).unwrap();

            let mut lck = stdout_ref.lock();

            writeln!(&mut lck, "{}\t{}\t{}", kmer, rvc, f).expect("Couldn't write output");
        }
    });*/
    let duration = start.elapsed();

    eprintln!(
        "Time elapsed creating hashmaps of all kmers in all sequences: {:?}\n",
        hash_duration
    );

    eprintln!("Time elapsed merging hashmaps: {:?}\n", uniq_duration);

    eprintln!("Time elapsed in runtime: {:?}\n", duration);

    Ok(())
}
