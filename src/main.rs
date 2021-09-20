use kmer_basic::Config;
use std::env;
use std::process;

fn main() {
    let config = Config::new(env::args()).unwrap_or_else(|err| {
        eprintln!("Problem parsing arguments: {}", err);
        process::exit(1);
    });

    eprintln!("\nSearching for kmers of length {}", config.kmer_len);
    eprintln!("... in file {}\n", config.filepath);

    if let Err(e) = kmer_basic::run(config) {
        eprintln!("Application error: {}", e);
        process::exit(1);
    }
}
