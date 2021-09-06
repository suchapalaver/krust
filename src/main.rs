use std::env;
use std::process;

use krust::Config;

fn main() {
   
    let config = Config::new(env::args()).unwrap_or_else(|err| {
        println!("Problem parsing arguments: {}", err);
        process::exit(1);
    });
    
    println!("\nSearching for kmers of length {}", config.kmer_len);
    println!("... in file {}\n", config.filepath);
    
    if let Err(e) = krust::run(config) {
        println!("Application error: {}", e);
        process::exit(1);
    }
}   
