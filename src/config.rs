//! CLI configuration and validation.
//!
//! This module handles parsing and validating command-line arguments.

use std::{error::Error, fs, path::PathBuf};

use colored::Colorize;

/// Configuration for k-mer counting.
#[derive(Debug)]
pub struct Config {
    /// The k-mer length (1-32).
    pub k: usize,
    /// Path to the input FASTA file.
    pub path: PathBuf,
}

impl Config {
    /// Creates a new configuration from string arguments.
    ///
    /// # Arguments
    ///
    /// * `k` - The k-mer length as a string (must be 1-32)
    /// * `path` - Path to a FASTA file (must exist)
    ///
    /// # Errors
    ///
    /// Returns an error if `k` is not a valid number in range 1-32,
    /// or if the file at `path` does not exist.
    pub fn new(k: &str, path: &str) -> Result<Config, Box<dyn Error>> {
        let k: usize = match k.parse::<usize>() {
            Ok(k) if k > 0 && k < 33 => k,
            Ok(_) => return Err("k-mer length needs to be larger than zero and, for krust currently, no more than 32".into()),
            Err(_) => return Err(format!("Issue with k-mer length argument \"{k}\"", k = k.bold()).into()),
        };

        let path = match fs::metadata(path) {
            Ok(_) => path.into(),
            Err(e) => {
                return Err(
                    format!("Issue with file path: {err}", err = e.to_string().bold()).into(),
                )
            }
        };

        Ok(Config { k, path })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    #[test]
    fn config_rejects_k_zero() {
        let file = NamedTempFile::new().unwrap();
        let result = Config::new("0", file.path().to_str().unwrap());
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("larger than zero"));
    }

    #[test]
    fn config_rejects_k_greater_than_32() {
        let file = NamedTempFile::new().unwrap();
        let result = Config::new("33", file.path().to_str().unwrap());
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("no more than 32"));
    }

    #[test]
    fn config_accepts_k_in_valid_range() {
        let file = NamedTempFile::new().unwrap();
        let path = file.path().to_str().unwrap();

        for k in 1..=32 {
            let result = Config::new(&k.to_string(), path);
            assert!(result.is_ok(), "k={k} should be valid");
            assert_eq!(result.unwrap().k, k);
        }
    }

    #[test]
    fn config_rejects_non_numeric_k() {
        let file = NamedTempFile::new().unwrap();
        let result = Config::new("abc", file.path().to_str().unwrap());
        assert!(result.is_err());
    }

    #[test]
    fn config_rejects_negative_k() {
        let file = NamedTempFile::new().unwrap();
        let result = Config::new("-5", file.path().to_str().unwrap());
        assert!(result.is_err());
    }

    #[test]
    fn config_rejects_invalid_path() {
        let result = Config::new("5", "/nonexistent/path/to/file.fa");
        assert!(result.is_err());
        assert!(result
            .unwrap_err()
            .to_string()
            .contains("Issue with file path"));
    }

    #[test]
    fn config_accepts_valid_file() {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(file, ">seq1\nACGT").unwrap();

        let result = Config::new("4", file.path().to_str().unwrap());
        assert!(result.is_ok());
        let config = result.unwrap();
        assert_eq!(config.k, 4);
        assert_eq!(config.path, file.path());
    }
}
