use std::{env, error::Error, path::PathBuf};

/// Parsing command line k-size and filepath arguments.
pub struct Config {
    pub k: usize,
    pub path: PathBuf,
}

impl Config {
    pub fn new(mut args: env::Args) -> Result<Config, Box<dyn Error>> {
        let k: usize = match args.nth(1) {
            Some(arg) => match arg.parse() {
                Ok(k) if k > 0 && k < 33 => k,
                Ok(_) => return Err("k-mer length needs to be larger than zero and, for `krust` in its current working form, no more than 32".into()),
                Err(_) => return Err(format!("issue with k-mer length argument: {}", arg).into()),
            },
            None => return Err("k-mer length input required".into()),
        };

        let path = match args.next() {
            Some(arg) => PathBuf::from(arg),
            None => return Err("filepath argument needed".into()),
        };

        Ok(Config { k, path })
    }
}
