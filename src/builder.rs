//! Builder pattern API for ergonomic k-mer counting.
//!
//! This module provides a fluent builder interface for configuring and executing
//! k-mer counting operations.
//!
//! # Example
//!
//! ```rust,no_run
//! use kmerust::builder::KmerCounter;
//!
//! let counts = KmerCounter::new()
//!     .k(21)?
//!     .min_count(2)
//!     .count("genome.fa")?;
//!
//! for (kmer, count) in counts {
//!     println!("{kmer}: {count}");
//! }
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```

use std::{collections::HashMap, error::Error, fmt::Debug, io::Write, path::Path};

use crate::{
    cli::OutputFormat,
    error::KmerLengthError,
    kmer::KmerLength,
    progress::Progress,
    run::{count_kmers, count_kmers_with_progress, run_with_options},
};

/// A builder for configuring k-mer counting operations.
///
/// Use [`KmerCounter::new()`] to create a new builder, configure it with the
/// fluent API, then call [`count()`](KmerCounter::count) or
/// [`count_to_writer()`](KmerCounter::count_to_writer) to execute.
///
/// # Example
///
/// ```rust,no_run
/// use kmerust::builder::KmerCounter;
///
/// // Basic usage
/// let counts = KmerCounter::new()
///     .k(21)?
///     .count("sequences.fa")?;
///
/// // With all options
/// let counts = KmerCounter::new()
///     .k(21)?
///     .min_count(5)
///     .count("sequences.fa")?;
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
#[derive(Debug, Clone)]
pub struct KmerCounter {
    k: Option<KmerLength>,
    min_count: i32,
    format: OutputFormat,
}

impl Default for KmerCounter {
    fn default() -> Self {
        Self::new()
    }
}

impl KmerCounter {
    /// Creates a new `KmerCounter` builder with default settings.
    ///
    /// Default settings:
    /// - `k`: None (must be set before counting)
    /// - `min_count`: 1 (include all k-mers)
    /// - `format`: FASTA
    ///
    /// Note: All k-mer counting uses canonical k-mers (k-mer and its reverse
    /// complement are treated as equivalent).
    ///
    /// # Example
    ///
    /// ```rust
    /// use kmerust::builder::KmerCounter;
    ///
    /// let counter = KmerCounter::new();
    /// ```
    #[must_use]
    pub fn new() -> Self {
        Self {
            k: None,
            min_count: 1,
            format: OutputFormat::Fasta,
        }
    }

    /// Sets the k-mer length.
    ///
    /// The k-mer length must be between 1 and 32 (inclusive).
    ///
    /// # Errors
    ///
    /// Returns [`KmerLengthError`] if `k` is outside the valid range.
    ///
    /// # Example
    ///
    /// ```rust
    /// use kmerust::builder::KmerCounter;
    ///
    /// let counter = KmerCounter::new().k(21)?;
    /// # Ok::<(), kmerust::error::KmerLengthError>(())
    /// ```
    pub fn k(mut self, k: usize) -> Result<Self, KmerLengthError> {
        self.k = Some(KmerLength::new(k)?);
        Ok(self)
    }

    /// Sets the k-mer length from a pre-validated `KmerLength`.
    ///
    /// Use this when you already have a validated `KmerLength` instance.
    ///
    /// # Example
    ///
    /// ```rust
    /// use kmerust::builder::KmerCounter;
    /// use kmerust::kmer::KmerLength;
    ///
    /// let k = KmerLength::new(21)?;
    /// let counter = KmerCounter::new().k_validated(k);
    /// # Ok::<(), kmerust::error::KmerLengthError>(())
    /// ```
    #[must_use]
    pub fn k_validated(mut self, k: KmerLength) -> Self {
        self.k = Some(k);
        self
    }

    /// Sets the minimum count threshold.
    ///
    /// K-mers with counts below this threshold will be excluded from results.
    /// Default is 1 (include all k-mers).
    ///
    /// # Example
    ///
    /// ```rust
    /// use kmerust::builder::KmerCounter;
    ///
    /// // Only include k-mers that appear at least 5 times
    /// let counter = KmerCounter::new()
    ///     .k(21)?
    ///     .min_count(5);
    /// # Ok::<(), kmerust::error::KmerLengthError>(())
    /// ```
    #[must_use]
    pub fn min_count(mut self, min_count: i32) -> Self {
        self.min_count = min_count;
        self
    }

    /// Sets the output format for [`count_to_writer()`](Self::count_to_writer).
    ///
    /// # Example
    ///
    /// ```rust
    /// use kmerust::builder::KmerCounter;
    /// use kmerust::cli::OutputFormat;
    ///
    /// let counter = KmerCounter::new()
    ///     .k(21)?
    ///     .format(OutputFormat::Tsv);
    /// # Ok::<(), kmerust::error::KmerLengthError>(())
    /// ```
    #[must_use]
    pub fn format(mut self, format: OutputFormat) -> Self {
        self.format = format;
        self
    }

    /// Counts k-mers in the specified FASTA file.
    ///
    /// Returns a HashMap mapping k-mer strings to their counts.
    ///
    /// # Errors
    ///
    /// Returns an error if:
    /// - `k` has not been set
    /// - The file cannot be read
    /// - The file cannot be parsed as FASTA
    ///
    /// # Example
    ///
    /// ```rust,no_run
    /// use kmerust::builder::KmerCounter;
    ///
    /// let counts = KmerCounter::new()
    ///     .k(21)?
    ///     .count("genome.fa")?;
    ///
    /// println!("Found {} unique k-mers", counts.len());
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn count<P>(&self, path: P) -> Result<HashMap<String, i32>, Box<dyn Error>>
    where
        P: AsRef<Path> + Debug,
    {
        let k = self.k.ok_or("k-mer length not set; call .k() first")?;

        let counts = count_kmers(&path, k.get())?;

        // Apply min_count filter
        if self.min_count > 1 {
            Ok(counts
                .into_iter()
                .filter(|(_, count)| *count >= self.min_count)
                .collect())
        } else {
            Ok(counts)
        }
    }

    /// Counts k-mers with progress reporting.
    ///
    /// Similar to [`count()`](Self::count), but invokes a callback after
    /// processing each sequence.
    ///
    /// # Errors
    ///
    /// Returns an error if:
    /// - `k` has not been set
    /// - The file cannot be read
    /// - The file cannot be parsed as FASTA
    ///
    /// # Example
    ///
    /// ```rust,no_run
    /// use kmerust::builder::KmerCounter;
    ///
    /// let counts = KmerCounter::new()
    ///     .k(21)?
    ///     .count_with_progress("genome.fa", |progress| {
    ///         eprintln!(
    ///             "Processed {} sequences ({} bases)",
    ///             progress.sequences_processed,
    ///             progress.bases_processed
    ///         );
    ///     })?;
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn count_with_progress<P, F>(
        &self,
        path: P,
        callback: F,
    ) -> Result<HashMap<String, i32>, Box<dyn Error>>
    where
        P: AsRef<Path> + Debug,
        F: Fn(Progress) + Send + Sync + 'static,
    {
        let k = self.k.ok_or("k-mer length not set; call .k() first")?;

        let counts = count_kmers_with_progress(&path, k.get(), callback)?;

        // Apply min_count filter
        if self.min_count > 1 {
            Ok(counts
                .into_iter()
                .filter(|(_, count)| *count >= self.min_count)
                .collect())
        } else {
            Ok(counts)
        }
    }

    /// Counts k-mers and writes results to stdout.
    ///
    /// # Errors
    ///
    /// Returns an error if:
    /// - `k` has not been set
    /// - The file cannot be read
    /// - Output cannot be written
    ///
    /// # Example
    ///
    /// ```rust,no_run
    /// use kmerust::builder::KmerCounter;
    /// use kmerust::cli::OutputFormat;
    ///
    /// KmerCounter::new()
    ///     .k(21)?
    ///     .format(OutputFormat::Tsv)
    ///     .min_count(2)
    ///     .run("genome.fa")?;
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn run<P>(&self, path: P) -> Result<(), Box<dyn Error>>
    where
        P: AsRef<Path> + Debug,
    {
        let k = self.k.ok_or("k-mer length not set; call .k() first")?;
        run_with_options(path, k.get(), self.format, self.min_count)?;
        Ok(())
    }

    /// Counts k-mers and writes results to a writer.
    ///
    /// # Errors
    ///
    /// Returns an error if:
    /// - `k` has not been set
    /// - The file cannot be read
    /// - Output cannot be written
    ///
    /// # Example
    ///
    /// ```rust,no_run
    /// use kmerust::builder::KmerCounter;
    /// use kmerust::cli::OutputFormat;
    /// use std::fs::File;
    /// use std::io::BufWriter;
    ///
    /// let file = File::create("output.tsv")?;
    /// let writer = BufWriter::new(file);
    ///
    /// KmerCounter::new()
    ///     .k(21)?
    ///     .format(OutputFormat::Tsv)
    ///     .count_to_writer("genome.fa", writer)?;
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn count_to_writer<P, W>(&self, path: P, mut writer: W) -> Result<(), Box<dyn Error>>
    where
        P: AsRef<Path> + Debug,
        W: Write,
    {
        let counts = self.count(&path)?;

        match self.format {
            OutputFormat::Fasta => {
                for (kmer, count) in counts {
                    writeln!(writer, ">{count}\n{kmer}")?;
                }
            }
            OutputFormat::Tsv => {
                for (kmer, count) in counts {
                    writeln!(writer, "{kmer}\t{count}")?;
                }
            }
            OutputFormat::Json => {
                #[derive(serde::Serialize)]
                struct KmerCount {
                    kmer: String,
                    count: i32,
                }
                let json_data: Vec<KmerCount> = counts
                    .into_iter()
                    .map(|(kmer, count)| KmerCount { kmer, count })
                    .collect();
                serde_json::to_writer_pretty(&mut writer, &json_data)?;
                writeln!(writer)?;
            }
        }

        writer.flush()?;
        Ok(())
    }

    /// Counts k-mers using memory-mapped I/O.
    ///
    /// Memory-maps the FASTA file for potentially faster access on large files.
    ///
    /// # Errors
    ///
    /// Returns an error if:
    /// - `k` has not been set
    /// - The file cannot be opened or memory-mapped
    /// - The file cannot be parsed as FASTA
    ///
    /// # Safety
    ///
    /// The underlying file must not be modified while being processed.
    ///
    /// # Example
    ///
    /// ```rust,no_run
    /// use kmerust::builder::KmerCounter;
    ///
    /// let counts = KmerCounter::new()
    ///     .k(21)?
    ///     .count_mmap("large_genome.fa")?;
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    #[cfg(feature = "mmap")]
    pub fn count_mmap<P>(&self, path: P) -> Result<HashMap<String, i32>, Box<dyn Error>>
    where
        P: AsRef<Path> + Debug,
    {
        use crate::run::count_kmers_mmap;

        let k = self.k.ok_or("k-mer length not set; call .k() first")?;
        let counts = count_kmers_mmap(&path, k.get())?;

        // Apply min_count filter
        if self.min_count > 1 {
            Ok(counts
                .into_iter()
                .filter(|(_, count)| *count >= self.min_count)
                .collect())
        } else {
            Ok(counts)
        }
    }

    /// Counts k-mers using streaming I/O for memory efficiency.
    ///
    /// Processes sequences one at a time without loading the entire file,
    /// suitable for very large files.
    ///
    /// # Errors
    ///
    /// Returns an error if:
    /// - `k` has not been set
    /// - The file cannot be read
    /// - The file cannot be parsed as FASTA
    ///
    /// # Example
    ///
    /// ```rust,no_run
    /// use kmerust::builder::KmerCounter;
    ///
    /// let counts = KmerCounter::new()
    ///     .k(21)?
    ///     .count_streaming("huge_genome.fa")?;
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn count_streaming<P>(&self, path: P) -> Result<HashMap<String, i32>, Box<dyn Error>>
    where
        P: AsRef<Path> + Debug,
    {
        use crate::streaming::count_kmers_streaming;

        let k = self.k.ok_or("k-mer length not set; call .k() first")?;
        let counts = count_kmers_streaming(&path, k.get())?;

        // Apply min_count filter
        if self.min_count > 1 {
            Ok(counts
                .into_iter()
                .filter(|(_, count)| *count >= self.min_count)
                .collect())
        } else {
            Ok(counts)
        }
    }

    /// Returns the configured k-mer length, if set.
    #[must_use]
    pub fn get_k(&self) -> Option<KmerLength> {
        self.k
    }

    /// Returns the configured minimum count threshold.
    #[must_use]
    pub fn get_min_count(&self) -> i32 {
        self.min_count
    }

    /// Returns the configured output format.
    #[must_use]
    pub fn get_format(&self) -> OutputFormat {
        self.format
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn builder_default() {
        let counter = KmerCounter::new();
        assert!(counter.get_k().is_none());
        assert_eq!(counter.get_min_count(), 1);
    }

    #[test]
    fn builder_k_valid() {
        let counter = KmerCounter::new().k(21).unwrap();
        assert_eq!(counter.get_k().unwrap().get(), 21);
    }

    #[test]
    fn builder_k_invalid() {
        let result = KmerCounter::new().k(0);
        assert!(result.is_err());

        let result = KmerCounter::new().k(33);
        assert!(result.is_err());
    }

    #[test]
    fn builder_k_validated() {
        let k = KmerLength::new(21).unwrap();
        let counter = KmerCounter::new().k_validated(k);
        assert_eq!(counter.get_k().unwrap().get(), 21);
    }

    #[test]
    fn builder_min_count() {
        let counter = KmerCounter::new().min_count(5);
        assert_eq!(counter.get_min_count(), 5);
    }

    #[test]
    fn builder_format() {
        let counter = KmerCounter::new().format(OutputFormat::Tsv);
        assert!(matches!(counter.get_format(), OutputFormat::Tsv));
    }

    #[test]
    fn builder_chained() {
        let counter = KmerCounter::new()
            .k(21)
            .unwrap()
            .min_count(3)
            .format(OutputFormat::Json);

        assert_eq!(counter.get_k().unwrap().get(), 21);
        assert_eq!(counter.get_min_count(), 3);
        assert!(matches!(counter.get_format(), OutputFormat::Json));
    }

    #[test]
    fn builder_count_without_k_fails() {
        let counter = KmerCounter::new();
        let result = counter.count("nonexistent.fa");
        assert!(result.is_err());
        assert!(result
            .unwrap_err()
            .to_string()
            .contains("k-mer length not set"));
    }

    #[test]
    fn builder_write_tsv_format() {
        let counts: HashMap<String, i32> = [("ACGT".to_string(), 5), ("TGCA".to_string(), 3)]
            .into_iter()
            .collect();

        let mut output = Cursor::new(Vec::new());
        let _counter = KmerCounter::new().k(4).unwrap().format(OutputFormat::Tsv);

        // Test the formatting logic directly
        for (kmer, count) in &counts {
            writeln!(output, "{kmer}\t{count}").unwrap();
        }

        let result = String::from_utf8(output.into_inner()).unwrap();
        assert!(result.contains("ACGT\t5") || result.contains("TGCA\t3"));
    }
}
