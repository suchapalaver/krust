use std::{error::Error, fs, path::PathBuf};

use colored::Colorize;

pub struct Config {
    pub k: usize,
    pub path: PathBuf,
}

impl Config {
    pub fn new(k: &str, path: &str) -> Result<Config, Box<dyn Error>> {
        let k: usize = match k.parse::<usize>() {
            Ok(k) if k > 0 && k < 33 => k,
            Ok(_) => return Err("k-mer length needs to be larger than zero and, for krust currently, no more than 32".into()),
            Err(_) => return Err(format!("Issue with k-mer length argument \"{}\"", k.bold()).into()),
        };

        let path = match fs::metadata(path) {
            Ok(_) => path.into(),
            Err(e) => return Err(format!("Issue with file path: {}", e.to_string().bold()).into()),
        };

        Ok(Config { k, path })
    }
}
