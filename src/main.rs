use std::{env, process};

fn main() {
    let config = krust::Config::new(env::args()).unwrap_or_else(|err| {
        eprintln!("Problem parsing arguments: {}", err);
        process::exit(1);
    });

    if let Err(e) = krust::canonicalize_kmers(config.filepath, config.kmer_len) {
        eprintln!("Application error: {}", e);
        drop(e);
        process::exit(1);
    }
}
