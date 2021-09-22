use bio::{alignment::sparse::hash_kmers, alphabets::dna::revcomp, io::fasta};
use fxhash::FxHashMap; // Stack O suggestions
use rayon::prelude::*;
use rayon::iter::ParallelBridge;
use std::{env, error::Error, fs::File, io::Write, str, time::Instant};
use dashmap::DashMap;

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

    let mut seq_map = FxHashMap::default();

    for (kmer, kmer_pos) in hash_kmers(result_data.seq(), k) { // the Vec<u32> that hash_kmers (rust-bio) is a list of indices of kmer's positions in seq. 
        seq_map.insert(kmer, kmer_pos.len());
    }
    seq_map
}

pub fn run(config: Config) -> Result<(), Box<dyn Error>> {
    let start = Instant::now();

    let filepath: String = config.filepath;

    let k: usize = config.kmer_len;

    let reader: fasta::Reader<std::io::BufReader<File>> =
        fasta::Reader::from_file(&filepath).unwrap();

    let fasta_records: Vec<Result<fasta::Record, std::io::Error>> = reader.records().collect();

    eprintln!("Number of records in fasta file: {}\n", fasta_records.len());
    
    let hash_vec: Vec<FxHashMap<&[u8], usize>> = fasta_records
        .par_iter() //  "Where you use par_iter(), instead of using map try using fold then reduce. With rayon, fold will let you merge the data into HashMaps in parallel. Then reduce will take those maps and let you merge them into a single one. If you want to avoid the Vec altogether, you should be able to call par_bridge directly on the records() result instead of calling collect (then call fold and reduce). par_bridge creates a parallel iterator from a regular iterator."
        .map(|result| hash_fasta_rec(result, k))
        .collect();

    let hash_duration = start.elapsed();

    eprintln!(
        "Time elapsed creating hashmaps of all kmers in all sequences: {:?}\n",
        hash_duration
    );
    
    eprintln!("number of hashmaps in vec: {}\n", hash_vec.len());

    // MERGING HASHMAPS

    // THIS WORKS
    let final_hash: DashMap<&[u8], usize> = DashMap::new();
    
    hash_vec.par_iter().for_each(|h| {
        for (kmer, freq) in h {
            if final_hash.contains_key(kmer) {
		*final_hash.get_mut(kmer).unwrap() += freq; // https://docs.rs/dashmap/4.0.2/dashmap/struct.DashMap.html
            } else {
                final_hash.insert(kmer, *freq);
            }
        }
    });
        
    // THIS WORKS
    /*
    let final_hash: FxHashMap<&[u8], usize> = hash_vec.par_iter()
    .fold(||FxHashMap::new(), |mut a: FxHashMap<&[u8], usize>, b| {
    a.extend(b.into_iter()); a}
)
    .reduce(||FxHashMap::new(),|mut a, b| {
    for (k, v) in b {
    if a.contains_key(k) {
    a.insert(k, v + a[k]);
} else {
    a.insert(k, v);
}
}
    a
}
);
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
    let uniq_duration = start.elapsed();

    eprintln!("Time elapsed merging hashmaps: {:?}\n", uniq_duration);
    
    // END OF MERGING

    let stdout_ref = &std::io::stdout();

    final_hash.into_iter().par_bridge().for_each(|(k, f)| {
        let kmer = str::from_utf8(k).unwrap();

        if kmer.contains("N") {
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
