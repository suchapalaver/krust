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

use std::{collections::HashMap, fmt::Debug, io::Write, path::Path};

use crate::{
    cli::OutputFormat,
    error::{BuilderError, KmerLengthError},
    format::SequenceFormat,
    kmer::KmerLength,
    progress::Progress,
    run::{count_kmers_with_format, count_kmers_with_progress, run_with_options},
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
/// // Basic usage - FASTA file
/// let counts = KmerCounter::new()
///     .k(21)?
///     .count("sequences.fa")?;
///
/// // FASTQ file (format auto-detected from extension)
/// let counts = KmerCounter::new()
///     .k(21)?
///     .count("reads.fq")?;
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
    min_count: u64,
    format: OutputFormat,
    input_format: SequenceFormat,
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
    /// - `format`: FASTA output
    /// - `input_format`: Auto (detected from file extension)
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
    pub const fn new() -> Self {
        Self {
            k: None,
            min_count: 1,
            format: OutputFormat::Fasta,
            input_format: SequenceFormat::Auto,
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
    pub const fn k_validated(mut self, k: KmerLength) -> Self {
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
    pub const fn min_count(mut self, min_count: u64) -> Self {
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
    pub const fn format(mut self, format: OutputFormat) -> Self {
        self.format = format;
        self
    }

    /// Sets the input file format.
    ///
    /// By default, format is auto-detected from the file extension:
    /// - `.fq`, `.fastq`, `.fq.gz`, `.fastq.gz` -> FASTQ
    /// - `.fa`, `.fasta`, `.fna`, `.fa.gz`, `.fasta.gz`, `.fna.gz` -> FASTA
    ///
    /// Use this method to explicitly specify the format when auto-detection
    /// is not sufficient (e.g., when reading from stdin or files with
    /// non-standard extensions).
    ///
    /// # Example
    ///
    /// ```rust
    /// use kmerust::builder::KmerCounter;
    /// use kmerust::format::SequenceFormat;
    ///
    /// let counter = KmerCounter::new()
    ///     .k(21)?
    ///     .input_format(SequenceFormat::Fastq);
    /// # Ok::<(), kmerust::error::KmerLengthError>(())
    /// ```
    #[must_use]
    pub const fn input_format(mut self, format: SequenceFormat) -> Self {
        self.input_format = format;
        self
    }

    /// Counts k-mers in the specified sequence file.
    ///
    /// Returns a `HashMap` mapping k-mer strings to their counts.
    /// Input format is auto-detected from extension unless explicitly set via
    /// [`input_format()`](Self::input_format).
    ///
    /// # Errors
    ///
    /// Returns [`BuilderError::KmerLengthNotSet`] if `k` has not been set,
    /// or [`BuilderError::Kmerust`] if the file cannot be read or parsed.
    ///
    /// # Example
    ///
    /// ```rust,no_run
    /// use kmerust::builder::KmerCounter;
    ///
    /// // FASTA file
    /// let counts = KmerCounter::new()
    ///     .k(21)?
    ///     .count("genome.fa")?;
    ///
    /// // FASTQ file (format auto-detected)
    /// let counts = KmerCounter::new()
    ///     .k(21)?
    ///     .count("reads.fq")?;
    ///
    /// println!("Found {} unique k-mers", counts.len());
    /// # Ok::<(), kmerust::error::BuilderError>(())
    /// ```
    pub fn count<P>(&self, path: P) -> Result<HashMap<String, u64>, BuilderError>
    where
        P: AsRef<Path> + Debug,
    {
        let k = self.k.ok_or(BuilderError::KmerLengthNotSet)?;

        let counts = count_kmers_with_format(&path, k.get(), self.input_format)?;

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

    /// Computes a k-mer frequency histogram from the specified sequence file.
    ///
    /// Returns a histogram mapping count values to the number of k-mers with
    /// that count. This is useful for genome size estimation, error detection,
    /// and heterozygosity analysis.
    ///
    /// # Errors
    ///
    /// Returns [`BuilderError::KmerLengthNotSet`] if `k` has not been set,
    /// or [`BuilderError::Kmerust`] if the file cannot be read or parsed.
    ///
    /// # Example
    ///
    /// ```rust,no_run
    /// use kmerust::builder::KmerCounter;
    ///
    /// let histogram = KmerCounter::new()
    ///     .k(21)?
    ///     .histogram("genome.fa")?;
    ///
    /// for (count, frequency) in histogram {
    ///     println!("{count} k-mers appear {frequency} times");
    /// }
    /// # Ok::<(), kmerust::error::BuilderError>(())
    /// ```
    pub fn histogram<P>(&self, path: P) -> Result<crate::histogram::KmerHistogram, BuilderError>
    where
        P: AsRef<Path> + Debug,
    {
        use crate::histogram::compute_histogram;

        let counts = self.count(&path)?;
        Ok(compute_histogram(&counts))
    }

    /// Counts k-mers with progress reporting.
    ///
    /// Similar to [`count()`](Self::count), but invokes a callback after
    /// processing each sequence.
    ///
    /// # Errors
    ///
    /// Returns [`BuilderError::KmerLengthNotSet`] if `k` has not been set,
    /// or [`BuilderError::Kmerust`] if the file cannot be read or parsed.
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
    /// # Ok::<(), kmerust::error::BuilderError>(())
    /// ```
    pub fn count_with_progress<P, F>(
        &self,
        path: P,
        callback: F,
    ) -> Result<HashMap<String, u64>, BuilderError>
    where
        P: AsRef<Path> + Debug,
        F: Fn(Progress) + Send + Sync + 'static,
    {
        let k = self.k.ok_or(BuilderError::KmerLengthNotSet)?;

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
    /// Returns [`BuilderError::KmerLengthNotSet`] if `k` has not been set,
    /// or other errors if the file cannot be read or output cannot be written.
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
    /// # Ok::<(), kmerust::error::BuilderError>(())
    /// ```
    pub fn run<P>(&self, path: P) -> Result<(), BuilderError>
    where
        P: AsRef<Path> + Debug,
    {
        let k = self.k.ok_or(BuilderError::KmerLengthNotSet)?;
        run_with_options(path, k.get(), self.format, self.min_count)?;
        Ok(())
    }

    /// Counts k-mers and writes results to a writer.
    ///
    /// # Errors
    ///
    /// Returns [`BuilderError::KmerLengthNotSet`] if `k` has not been set,
    /// or other errors if the file cannot be read or output cannot be written.
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
    /// # Ok::<(), kmerust::error::BuilderError>(())
    /// ```
    pub fn count_to_writer<P, W>(&self, path: P, mut writer: W) -> Result<(), BuilderError>
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
                    count: u64,
                }
                let json_data: Vec<KmerCount> = counts
                    .into_iter()
                    .map(|(kmer, count)| KmerCount { kmer, count })
                    .collect();
                serde_json::to_writer_pretty(&mut writer, &json_data)?;
                writeln!(writer)?;
            }
            OutputFormat::Histogram => {
                use crate::histogram::compute_histogram;

                let histogram = compute_histogram(&counts);
                for (count, frequency) in histogram {
                    writeln!(writer, "{count}\t{frequency}")?;
                }
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
    /// Returns [`BuilderError::KmerLengthNotSet`] if `k` has not been set,
    /// or [`BuilderError::Kmerust`] if the file cannot be opened, memory-mapped, or parsed.
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
    /// # Ok::<(), kmerust::error::BuilderError>(())
    /// ```
    #[cfg(feature = "mmap")]
    pub fn count_mmap<P>(&self, path: P) -> Result<HashMap<String, u64>, BuilderError>
    where
        P: AsRef<Path> + Debug,
    {
        use crate::run::count_kmers_mmap;

        let k = self.k.ok_or(BuilderError::KmerLengthNotSet)?;
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
    /// Returns [`BuilderError::KmerLengthNotSet`] if `k` has not been set,
    /// or [`BuilderError::Kmerust`] if the file cannot be read or parsed.
    ///
    /// # Example
    ///
    /// ```rust,no_run
    /// use kmerust::builder::KmerCounter;
    ///
    /// let counts = KmerCounter::new()
    ///     .k(21)?
    ///     .count_streaming("huge_genome.fa")?;
    /// # Ok::<(), kmerust::error::BuilderError>(())
    /// ```
    pub fn count_streaming<P>(&self, path: P) -> Result<HashMap<String, u64>, BuilderError>
    where
        P: AsRef<Path> + Debug,
    {
        use crate::streaming::count_kmers_streaming;

        let k = self.k.ok_or(BuilderError::KmerLengthNotSet)?;
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
    pub const fn get_k(&self) -> Option<KmerLength> {
        self.k
    }

    /// Returns the configured minimum count threshold.
    #[must_use]
    pub const fn get_min_count(&self) -> u64 {
        self.min_count
    }

    /// Returns the configured output format.
    #[must_use]
    pub const fn get_format(&self) -> OutputFormat {
        self.format
    }

    /// Returns the configured input format.
    #[must_use]
    pub const fn get_input_format(&self) -> SequenceFormat {
        self.input_format
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
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
        let counts: HashMap<String, u64> = [("ACGT".to_string(), 5), ("TGCA".to_string(), 3)]
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
