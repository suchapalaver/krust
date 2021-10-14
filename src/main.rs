use krust::{canonicalize_kmers, Config};
use std::{
    env,
    io::{BufWriter, Write},
    process, str,
};

fn main() {
    let config = Config::new(env::args()).unwrap_or_else(|err| {
        eprintln!("Problem parsing arguments: {}", err);
        process::exit(1);
    });

    assert!(
        config.kmer_len > 0,
        "The requested k-mer length caused a problem. k = {}",
        config.kmer_len
    );

    eprintln!(
        "Searching for k-mers of length {}\n... in file {}\n",
        config.kmer_len, config.filepath
    );

    match canonicalize_kmers(config.filepath, config.kmer_len) {
        Err(e) => {
            eprintln!("Application error: {}", e);
            process::exit(1);
        }
        Ok(fasta_hash) => {
            let mut buf = BufWriter::new(std::io::stdout());

            fasta_hash.into_iter().for_each(|(kmer, count)| {
                writeln!(buf, ">{}\n{}", count, str::from_utf8(&kmer).unwrap())
                    .expect("Unable to write output");
            });

            buf.flush().unwrap();
        }
    }
}
