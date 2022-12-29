use std::{error::Error, path::PathBuf, fs};

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
            Err(_) => return Err(format!("Issue with k-mer length argument \"{}\"", k).into()),
        };

        let path = match fs::metadata(path) {
            Ok(_) => path.into(),
            Err(e) => return Err(format!("Issue with file path: {}", e).into()),
        };

        let reader = match reader {
            reader if matches!(reader, "needletail") => true,
            reader if matches!(reader, "rust-bio") => true,
            _ => return Err(format!("Invalid reader argument: \"{}\"", reader).into()),
        };

        Ok(Config { k, path, reader })
    }
}
