use std::env;
use std::error::Error;

/// Parsing command line k-size and filepath arguments.
pub struct Config {
    pub kmer_len: usize,
    pub filepath: String,
}

impl Config {
    pub fn new(mut args: env::Args) -> Result<Config, Box<dyn Error>> {
        let kmer_len: usize = match args.nth(1) {
            Some(arg) => match arg.parse() {
                Ok(kmer_len) if kmer_len > 0 && kmer_len < 33 => kmer_len,
                Ok(_) => return Err("k-mer length needs to be larger than zero and, for `krust` in its current working form, no more than 32".into()),
                Err(_) => return Err(format!("issue with k-mer length argument: {}", arg).into()),
            },
            None => return Err("k-mer length input required".into()),
        };

        let filepath = match args.next() {
            Some(arg) => arg,
            None => return Err("filepath argument needed".into()),
        };

        Ok(Config { kmer_len, filepath })
    }
}
