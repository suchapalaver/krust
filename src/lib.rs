use bio::{alignment::sparse::hash_kmers, alphabets::dna::revcomp, io::fasta};
use itertools::Itertools;
use rayon::prelude::*;
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
    
    let hash_vec: Vec<HashMap<&[u8], usize>> = fasta_records
        .par_iter() //  Where you use par_iter(), instead of using map try using fold then reduce. With rayon, fold will let you merge the data into HashMaps in parallel. Then reduce will take those maps and let you merge them into a single one. If you want to avoid the Vec altogether, you should be able to call par_bridge directly on the records() result instead of calling collect (then call fold and reduce). par_bridge creates a parallel iterator from a regular iterator.
        .map(|result| hash_fasta_rec(result, k))
        .collect();
    
    let hash_duration = start.elapsed();

    eprintln!(
        "Time elapsed creating hashmaps of all kmers in all sequences: {:?}\n",
        hash_duration
    );
    // merging hashmaps
    //eprintln!("length of hash_vec now: {}", hash_vec.len());
    
//  Idea was to find the longest map and somehow make that the base for comparison but scrapping for now...
    /*
    let mut hash_len_vec = HashSet::new(); // create set of number of kmers 
    
    for h in &hash_vec {
	hash_len_vec.insert(h.len());
    }
    //eprintln!("hashmap lengths: {:?}", hash_len_vec);

    let longest_len = hash_len_vec.iter().max().unwrap();
    
    let i = &hash_vec.iter().position(|h| h.len() == *longest_len).unwrap();

    let mut final_hash = hash_vec.remove(*i);

    //eprintln!("this is the hash we're basing off: {:?}", final_hash);

    //eprintln!("length of hash_vec post removal: {}", hash_vec.len());
     */
    
/*
    let pair_hash = hash_vec.clone();
    
    let mut pairing_hash: HashMap<&[u8], usize> = HashMap::new();
    
    for pair in &pair_hash.into_iter().chunks(2) {
	for p in pair {
	    eprintln!("p.len(): {}", p.len());
	    if pairing_hash.is_empty() {
		let pairing_hash = p;
		continue;
	    } else {
		for (kmer, freq) in p {
		    if pairing_hash.contains_key(kmer) {
			pairing_hash.insert(kmer, freq + pairing_hash[kmer]);
		    } else {
			pairing_hash.insert(kmer, freq);
		    }
		}
	    }
	}
    }
    eprintln!("{}", pairing_hash.len());
     */
    /*
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
    

    eprintln!("number of hashmaps in vec: {}", hash_vec.len());

    //let mut final_hash = HashMap::new();

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
    
    //let merged: HashMap<K, V> = HashMap::new();
	
	/*
	    for (kmer, freq) in pair {
		if merged.contains_key(kmer) {
                    merged.insert(kmer, freq + merged[kmer]);
		} else {
                    merged.insert(kmer, freq);
		}
            }
	}
    };
    */
    
    let uniq_duration = start.elapsed();

    eprintln!("Time elapsed merging hashmaps: {:?}\n", uniq_duration);

    let stdout_ref = &std::io::stdout();

    final_hash.par_iter().for_each(|(k, f)| {
        let kmer = str::from_utf8(k).unwrap();

        if kmer.contains("N") {
        } else {
            let rvc = revcomp(*k);

            let rvc = str::from_utf8(&rvc).unwrap();

            let mut lck = stdout_ref.lock();

            writeln!(&mut lck, "{}\t{}\t{}", kmer, rvc, f).expect("Couldn't write output");
        }
    });
    let duration = start.elapsed();

    eprintln!("Time elapsed in runtime: {:?}\n", duration);

    Ok(())
}
