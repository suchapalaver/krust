use std::{error::Error, path::PathBuf};

pub struct Config {
    pub k: usize,
    pub path: PathBuf,
    pub reader: bool,
}

impl Config {
    pub fn new(k: &str, path: &str, reader: &str) -> Result<Config, Box<dyn Error>> {
        let k: usize = match k.parse::<usize>() {
            Ok(k) if k > 0 && k < 33 => k,
            Ok(_) => return Err("k-mer length needs to be larger than zero and, for krust currently, no more than 32".into()),
            Err(_) => return Err(format!("issue with k-mer length argument: {}", k).into()),
        };

        let path = path.into();

        let reader = matches!(reader, "needletail");

        Ok(Config { k, path, reader })
    }
}
