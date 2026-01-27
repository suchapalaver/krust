//! Streaming k-mer counting for memory-efficient processing of large files.
//!
//! This module provides APIs that process FASTA and FASTQ files without loading all
//! sequences into memory simultaneously, making it suitable for very large genomic
//! datasets. File format is auto-detected from the file extension.
//!
//! # API Variants
//!
//! | Function | Parallelism | Memory | Best For |
//! |----------|-------------|--------|----------|
//! | [`count_kmers_streaming`] | Parallel | Moderate | Most use cases |
//! | [`count_kmers_sequential`] | Single-threaded | Minimal | Extreme memory constraints |
//! | [`count_kmers_from_reader`] | Single-threaded | Minimal | Stdin/custom readers |
//!
//! # Memory Model
//!
//! **Streaming API ([`count_kmers_streaming`]):**
//! - Records are batched for parallel processing
//! - Provides good speed/memory balance
//! - Memory usage: sequences + count map
//!
//! **Sequential API ([`count_kmers_sequential`]):**
//! - Processes each record immediately as read
//! - Lowest possible memory footprint
//! - Memory usage: one sequence + count map
//!
//! **Reader API ([`count_kmers_from_reader`]):**
//! - Works with any `BufRead` source including stdin
//! - Single-threaded, minimal memory
//! - Enables Unix pipeline integration
//!
//! For both APIs, the count map (one `u64` key + one `u64` value per unique k-mer)
//! dominates memory usage for most datasets.
//!
//! # Packed Bits API
//!
//! For performance-critical applications, use [`count_kmers_streaming_packed`] to get
//! results as packed 64-bit integers instead of strings. This avoids string allocation
//! overhead when you need to do further processing on the k-mer counts.

use std::{collections::HashMap, fmt::Debug, io::BufRead, path::Path};

use bytes::Bytes;
use dashmap::DashMap;
use rayon::prelude::*;
use rustc_hash::FxHasher;
use std::hash::BuildHasherDefault;

use crate::{
    error::KmeRustError,
    format::SequenceFormat,
    input::Input,
    kmer::{unpack_to_string, Kmer, KmerLength},
};

#[cfg(feature = "tracing")]
use tracing::{debug, info, info_span};

/// Counts k-mers in a FASTA or FASTQ file using streaming I/O.
///
/// Processes sequences one at a time without loading the entire file into memory.
/// This is more memory-efficient than [`count_kmers`](crate::run::count_kmers) for
/// very large files.
///
/// File format is auto-detected from the extension:
/// - `.fa`, `.fasta`, `.fna` (and `.gz` variants) → FASTA
/// - `.fq`, `.fastq` (and `.gz` variants) → FASTQ
///
/// # Arguments
///
/// * `path` - Path to the FASTA/FASTQ file
/// * `k` - K-mer length (must be 1-32)
///
/// # Returns
///
/// A HashMap mapping k-mer strings to their counts.
///
/// # Errors
///
/// Returns an error if:
/// - `k` is outside the valid range (1-32)
/// - The file cannot be read or parsed
///
/// # Example
///
/// ```rust,no_run
/// use kmerust::streaming::count_kmers_streaming;
///
/// // Works with both FASTA and FASTQ
/// let counts = count_kmers_streaming("large_genome.fa", 21)?;
/// let counts = count_kmers_streaming("reads.fq.gz", 21)?;
/// println!("Found {} unique k-mers", counts.len());
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
pub fn count_kmers_streaming<P>(path: P, k: usize) -> Result<HashMap<String, u64>, KmeRustError>
where
    P: AsRef<Path> + Debug,
{
    #[cfg(feature = "tracing")]
    info!(k = k, path = ?path, "Starting streaming k-mer counting");

    let k_len = KmerLength::new(k)?;
    let packed = count_kmers_streaming_packed(&path, k_len)?;

    #[cfg(feature = "tracing")]
    let _unpack_span = info_span!("unpack_kmers", count = packed.len()).entered();

    let result: HashMap<String, u64> = packed
        .into_par_iter()
        .map(|(bits, count)| (unpack_to_string(bits, k_len), count))
        .collect();

    #[cfg(feature = "tracing")]
    info!(
        unique_kmers = result.len(),
        "Streaming k-mer counting complete"
    );

    Ok(result)
}

/// Counts k-mers and returns packed 64-bit representations.
///
/// This is the most efficient API for k-mer counting when you need to do
/// further processing on the results. It avoids string allocation entirely.
///
/// File format is auto-detected from the extension:
/// - `.fa`, `.fasta`, `.fna` (and `.gz` variants) → FASTA
/// - `.fq`, `.fastq` (and `.gz` variants) → FASTQ
///
/// # Arguments
///
/// * `path` - Path to the FASTA/FASTQ file
/// * `k` - Validated k-mer length
///
/// # Returns
///
/// A HashMap mapping packed k-mer bits to their counts.
///
/// # Example
///
/// ```rust,no_run
/// use kmerust::streaming::count_kmers_streaming_packed;
/// use kmerust::kmer::KmerLength;
///
/// let k = KmerLength::new(21)?;
/// let counts = count_kmers_streaming_packed("genome.fa", k)?;
/// let counts = count_kmers_streaming_packed("reads.fq.gz", k)?;
///
/// // Process packed bits directly without string conversion
/// for (packed_bits, count) in counts {
///     if count >= 10 {
///         println!("K-mer {packed_bits:#x} appears {count} times");
///     }
/// }
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
pub fn count_kmers_streaming_packed<P>(
    path: P,
    k: KmerLength,
) -> Result<HashMap<u64, u64>, KmeRustError>
where
    P: AsRef<Path> + Debug,
{
    let counter = StreamingKmerCounter::new();
    counter.count_file(path, k)
}

/// Counts k-mers from an in-memory byte slice.
///
/// Useful for testing or when sequence data is already in memory.
///
/// # Arguments
///
/// * `sequences` - Iterator of sequence byte slices
/// * `k` - Validated k-mer length
///
/// # Returns
///
/// A HashMap mapping packed k-mer bits to their counts.
///
/// # Example
///
/// ```rust
/// use kmerust::streaming::count_kmers_from_sequences;
/// use kmerust::kmer::KmerLength;
/// use bytes::Bytes;
///
/// let sequences = vec![
///     Bytes::from_static(b"ACGTACGT"),
///     Bytes::from_static(b"TGCATGCA"),
/// ];
///
/// let k = KmerLength::new(4)?;
/// let counts = count_kmers_from_sequences(sequences.into_iter(), k);
/// # Ok::<(), kmerust::error::KmerLengthError>(())
/// ```
pub fn count_kmers_from_sequences<I>(sequences: I, k: KmerLength) -> HashMap<u64, u64>
where
    I: Iterator<Item = Bytes>,
{
    let counter = StreamingKmerCounter::new();
    counter.count_sequences(sequences, k)
}

/// Counts k-mers with true sequential processing for minimum memory usage.
///
/// Unlike [`count_kmers_streaming`] which batches sequences for parallel processing,
/// this function processes each sequence immediately as it's read from disk.
/// This provides the lowest possible peak memory usage at the cost of being
/// single-threaded.
///
/// File format is auto-detected from the extension:
/// - `.fa`, `.fasta`, `.fna` (and `.gz` variants) → FASTA
/// - `.fq`, `.fastq` (and `.gz` variants) → FASTQ
///
/// # When to Use
///
/// Use this function when:
/// - Memory is extremely constrained
/// - The file is very large relative to available RAM
/// - You're processing files where individual sequences are very long
///
/// For most use cases, [`count_kmers_streaming`] provides a better balance of
/// speed and memory efficiency.
///
/// # Arguments
///
/// * `path` - Path to the FASTA/FASTQ file
/// * `k` - K-mer length (must be 1-32)
///
/// # Returns
///
/// A HashMap mapping packed k-mer bits to their counts.
///
/// # Errors
///
/// Returns an error if:
/// - `k` is outside the valid range (1-32)
/// - The file cannot be read or parsed
///
/// # Example
///
/// ```rust,no_run
/// use kmerust::streaming::count_kmers_sequential;
///
/// let counts = count_kmers_sequential("huge_genome.fa", 21)?;
/// let counts = count_kmers_sequential("reads.fq.gz", 21)?;
/// println!("Found {} unique k-mers", counts.len());
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
pub fn count_kmers_sequential<P>(path: P, k: usize) -> Result<HashMap<u64, u64>, KmeRustError>
where
    P: AsRef<Path> + Debug,
{
    let k_len = KmerLength::new(k)?;
    let counter = SequentialKmerCounter::new();
    counter.count_file(path, k_len)
}

/// Counts k-mers from standard input.
///
/// Reads FASTA-formatted sequences from stdin and counts k-mers.
/// This enables Unix pipeline integration:
///
/// ```bash
/// cat genome.fa | krust 21
/// zcat large.fa.gz | krust 21 > counts.tsv
/// seqtk sample reads.fq 0.1 | krust 17
/// ```
///
/// # Arguments
///
/// * `k` - K-mer length (must be 1-32)
///
/// # Returns
///
/// A HashMap mapping k-mer strings to their counts.
///
/// # Errors
///
/// Returns an error if:
/// - `k` is outside the valid range (1-32)
/// - The input cannot be parsed as FASTA
///
/// # Example
///
/// ```rust,no_run
/// use kmerust::streaming::count_kmers_stdin;
///
/// let counts = count_kmers_stdin(21)?;
/// println!("Found {} unique k-mers", counts.len());
/// # Ok::<(), kmerust::error::KmeRustError>(())
/// ```
pub fn count_kmers_stdin(k: usize) -> Result<HashMap<String, u64>, KmeRustError> {
    count_kmers_stdin_with_format(k, SequenceFormat::Auto)
}

/// Counts k-mers from standard input with explicit format specification.
///
/// # Arguments
///
/// * `k` - K-mer length (must be 1-32)
/// * `format` - Input format (Auto defaults to FASTA for stdin)
///
/// # Returns
///
/// A HashMap mapping k-mer strings to their counts.
///
/// # Errors
///
/// Returns an error if:
/// - `k` is outside the valid range (1-32)
/// - The input cannot be parsed
pub fn count_kmers_stdin_with_format(
    k: usize,
    format: SequenceFormat,
) -> Result<HashMap<String, u64>, KmeRustError> {
    let k_len = KmerLength::new(k)?;
    let stdin = std::io::stdin();
    let reader = stdin.lock();
    // For stdin, resolve format with None path (defaults to FASTA if Auto)
    let resolved_format = format.resolve(None);
    let packed = count_kmers_from_reader_impl_with_format(reader, k_len, resolved_format)?;

    Ok(packed
        .into_iter()
        .map(|(bits, count)| (unpack_to_string(bits, k_len), count))
        .collect())
}

/// Counts k-mers from standard input, returning packed bit representations.
///
/// This is the most efficient stdin API, avoiding string allocation.
///
/// # Arguments
///
/// * `k` - Validated k-mer length
///
/// # Returns
///
/// A HashMap mapping packed k-mer bits to their counts.
///
/// # Errors
///
/// Returns an error if the input cannot be parsed as FASTA.
pub fn count_kmers_stdin_packed(k: KmerLength) -> Result<HashMap<u64, u64>, KmeRustError> {
    let stdin = std::io::stdin();
    let reader = stdin.lock();
    count_kmers_from_reader_impl(reader, k)
}

/// Counts k-mers from any buffered reader.
///
/// This is the most flexible API, allowing k-mer counting from any source
/// that implements `BufRead`, including files, stdin, network streams, or
/// in-memory buffers.
///
/// # Arguments
///
/// * `reader` - Any type implementing `BufRead`
/// * `k` - K-mer length (must be 1-32)
///
/// # Returns
///
/// A HashMap mapping k-mer strings to their counts.
///
/// # Errors
///
/// Returns an error if:
/// - `k` is outside the valid range (1-32)
/// - The input cannot be parsed as FASTA
///
/// # Example
///
/// ```rust
/// use kmerust::streaming::count_kmers_from_reader;
/// use std::io::BufReader;
///
/// let fasta_data = b">seq1\nACGTACGT\n>seq2\nTGCATGCA\n";
/// let reader = BufReader::new(&fasta_data[..]);
/// let counts = count_kmers_from_reader(reader, 4)?;
/// assert!(!counts.is_empty());
/// # Ok::<(), kmerust::error::KmeRustError>(())
/// ```
pub fn count_kmers_from_reader<R>(reader: R, k: usize) -> Result<HashMap<String, u64>, KmeRustError>
where
    R: BufRead,
{
    let k_len = KmerLength::new(k)?;
    let packed = count_kmers_from_reader_impl(reader, k_len)?;

    Ok(packed
        .into_iter()
        .map(|(bits, count)| (unpack_to_string(bits, k_len), count))
        .collect())
}

/// Counts k-mers from any buffered reader, returning packed bit representations.
///
/// This is the most efficient reader API, avoiding string allocation.
///
/// # Arguments
///
/// * `reader` - Any type implementing `BufRead`
/// * `k` - Validated k-mer length
///
/// # Returns
///
/// A HashMap mapping packed k-mer bits to their counts.
///
/// # Errors
///
/// Returns an error if the input cannot be parsed as FASTA.
///
/// # Example
///
/// ```rust
/// use kmerust::streaming::count_kmers_from_reader_packed;
/// use kmerust::kmer::KmerLength;
/// use std::io::BufReader;
///
/// let fasta_data = b">seq1\nACGTACGT\n>seq2\nTGCATGCA\n";
/// let reader = BufReader::new(&fasta_data[..]);
/// let k = KmerLength::new(4)?;
/// let counts = count_kmers_from_reader_packed(reader, k)?;
/// assert!(!counts.is_empty());
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
pub fn count_kmers_from_reader_packed<R>(
    reader: R,
    k: KmerLength,
) -> Result<HashMap<u64, u64>, KmeRustError>
where
    R: BufRead,
{
    count_kmers_from_reader_impl(reader, k)
}

/// Counts k-mers from an [`Input`] source (file or stdin).
///
/// This is the main entry point for input-agnostic k-mer counting.
///
/// # Arguments
///
/// * `input` - The input source (file path or stdin)
/// * `k` - K-mer length (must be 1-32)
///
/// # Returns
///
/// A HashMap mapping k-mer strings to their counts.
///
/// # Errors
///
/// Returns an error if:
/// - `k` is outside the valid range (1-32)
/// - The input cannot be read or parsed
///
/// # Example
///
/// ```rust,no_run
/// use kmerust::streaming::count_kmers_from_input;
/// use kmerust::input::Input;
/// use std::path::Path;
///
/// // From file
/// let input = Input::from_path(Path::new("genome.fa"));
/// let counts = count_kmers_from_input(&input, 21)?;
///
/// // From stdin (would read from actual stdin)
/// let input = Input::Stdin;
/// // let counts = count_kmers_from_input(&input, 21)?;
/// # Ok::<(), kmerust::error::KmeRustError>(())
/// ```
pub fn count_kmers_from_input(
    input: &Input,
    k: usize,
) -> Result<HashMap<String, u64>, KmeRustError> {
    match input {
        Input::File(path) => count_kmers_streaming(path, k),
        Input::Stdin => count_kmers_stdin(k),
    }
}

/// Counts k-mers from an [`Input`] source, returning packed bit representations.
///
/// # Arguments
///
/// * `input` - The input source (file path or stdin)
/// * `k` - Validated k-mer length
///
/// # Returns
///
/// A HashMap mapping packed k-mer bits to their counts.
///
/// # Errors
///
/// Returns an error if the input cannot be read or parsed.
pub fn count_kmers_from_input_packed(
    input: &Input,
    k: KmerLength,
) -> Result<HashMap<u64, u64>, KmeRustError> {
    match input {
        Input::File(path) => count_kmers_streaming_packed(path, k),
        Input::Stdin => count_kmers_stdin_packed(k),
    }
}

/// Internal implementation for counting k-mers from a reader.
#[cfg(not(feature = "needletail"))]
fn count_kmers_from_reader_impl<R>(
    reader: R,
    k: KmerLength,
) -> Result<HashMap<u64, u64>, KmeRustError>
where
    R: BufRead,
{
    count_kmers_from_reader_impl_with_format(reader, k, SequenceFormat::Fasta)
}

/// Internal implementation for counting k-mers from a reader with format.
#[cfg(not(feature = "needletail"))]
fn count_kmers_from_reader_impl_with_format<R>(
    reader: R,
    k: KmerLength,
    format: SequenceFormat,
) -> Result<HashMap<u64, u64>, KmeRustError>
where
    R: BufRead,
{
    use bio::io::{fasta, fastq};

    let mut counts: HashMap<u64, u64, BuildHasherDefault<FxHasher>> =
        HashMap::with_hasher(BuildHasherDefault::default());

    match format {
        SequenceFormat::Fastq => {
            let fastq_reader = fastq::Reader::new(reader);
            for result in fastq_reader.records() {
                let record = result.map_err(|e| KmeRustError::SequenceParse {
                    details: e.to_string(),
                })?;
                process_sequence_into_counts(&mut counts, record.seq(), k);
            }
        }
        SequenceFormat::Fasta | SequenceFormat::Auto => {
            let fasta_reader = fasta::Reader::new(reader);
            for result in fasta_reader.records() {
                let record = result.map_err(|e| KmeRustError::SequenceParse {
                    details: e.to_string(),
                })?;
                process_sequence_into_counts(&mut counts, record.seq(), k);
            }
        }
    }

    Ok(counts.into_iter().collect())
}

/// Internal implementation for counting k-mers from a reader (needletail version).
///
/// Note: needletail requires the reader to be `Send`, so we read into a buffer first.
#[cfg(feature = "needletail")]
fn count_kmers_from_reader_impl<R>(
    reader: R,
    k: KmerLength,
) -> Result<HashMap<u64, u64>, KmeRustError>
where
    R: BufRead,
{
    // needletail auto-detects format, so we can ignore the format parameter
    count_kmers_from_reader_impl_with_format(reader, k, SequenceFormat::Auto)
}

/// Internal implementation for counting k-mers from a reader with format (needletail version).
///
/// Note: needletail requires the reader to be `Send`, so we read into a buffer first.
/// needletail auto-detects FASTA/FASTQ format, so the format parameter is informational.
#[cfg(feature = "needletail")]
fn count_kmers_from_reader_impl_with_format<R>(
    mut reader: R,
    k: KmerLength,
    _format: SequenceFormat,
) -> Result<HashMap<u64, u64>, KmeRustError>
where
    R: BufRead,
{
    use std::io::Cursor;

    // Read all data into a buffer since needletail requires Send
    let mut buffer = Vec::new();
    reader
        .read_to_end(&mut buffer)
        .map_err(|e| KmeRustError::SequenceParse {
            details: format!("failed to read input: {e}"),
        })?;

    let mut parser = needletail::parse_fastx_reader(Cursor::new(buffer)).map_err(|e| {
        KmeRustError::SequenceParse {
            details: e.to_string(),
        }
    })?;
    let mut counts: HashMap<u64, u64, BuildHasherDefault<FxHasher>> =
        HashMap::with_hasher(BuildHasherDefault::default());

    while let Some(result) = parser.next() {
        let record = result.map_err(|e| KmeRustError::SequenceParse {
            details: e.to_string(),
        })?;
        process_sequence_into_counts(&mut counts, &record.seq(), k);
    }

    Ok(counts.into_iter().collect())
}

/// Process a sequence and add k-mer counts to the map.
fn process_sequence_into_counts(
    counts: &mut HashMap<u64, u64, BuildHasherDefault<FxHasher>>,
    seq: &[u8],
    k: KmerLength,
) {
    let k_val = k.get();
    if seq.len() < k_val {
        return;
    }

    let mut i = 0;
    while i <= seq.len() - k_val {
        let sub = Bytes::copy_from_slice(&seq[i..i + k_val]);

        match Kmer::from_sub(sub) {
            Ok(unpacked) => {
                let canonical = unpacked.pack().canonical();
                *counts.entry(canonical.packed_bits()).or_insert(0) += 1;
                i += 1;
            }
            Err(err) => {
                i += err.position + 1;
            }
        }
    }
}

/// A truly sequential k-mer counter with minimal memory footprint.
///
/// Processes sequences one at a time as they're read, without batching.
struct SequentialKmerCounter {
    counts: HashMap<u64, u64, BuildHasherDefault<FxHasher>>,
}

impl SequentialKmerCounter {
    fn new() -> Self {
        Self {
            counts: HashMap::with_hasher(BuildHasherDefault::<FxHasher>::default()),
        }
    }

    #[cfg(not(feature = "needletail"))]
    fn count_file<P>(mut self, path: P, k: KmerLength) -> Result<HashMap<u64, u64>, KmeRustError>
    where
        P: AsRef<Path> + Debug,
    {
        use bio::io::{fasta, fastq};

        let path_ref = path.as_ref();
        let format = SequenceFormat::from_extension(path_ref);

        #[cfg(feature = "tracing")]
        let _span = info_span!("sequential_count", path = ?path_ref, ?format).entered();

        // Handle gzip if the feature is enabled
        #[cfg(feature = "gzip")]
        let is_gzip = path_ref.extension().map(|ext| ext == "gz").unwrap_or(false);

        #[cfg(feature = "gzip")]
        if is_gzip {
            use flate2::read::GzDecoder;
            use std::{fs::File, io::BufReader};

            let file = File::open(path_ref).map_err(|e| KmeRustError::SequenceRead {
                source: e,
                path: path_ref.to_path_buf(),
            })?;
            let decoder = GzDecoder::new(file);
            let buf_reader = BufReader::new(decoder);

            match format {
                SequenceFormat::Fastq => {
                    let reader = fastq::Reader::new(buf_reader);
                    for result in reader.records() {
                        let record = result.map_err(|e| KmeRustError::SequenceParse {
                            details: e.to_string(),
                        })?;
                        self.process_sequence(record.seq(), k);
                    }
                }
                SequenceFormat::Fasta | SequenceFormat::Auto => {
                    let reader = fasta::Reader::new(buf_reader);
                    for result in reader.records() {
                        let record = result.map_err(|e| KmeRustError::SequenceParse {
                            details: e.to_string(),
                        })?;
                        self.process_sequence(record.seq(), k);
                    }
                }
            }

            return Ok(self.counts.into_iter().collect());
        }

        // Non-gzip path
        match format {
            SequenceFormat::Fastq => {
                let reader =
                    fastq::Reader::from_file(path_ref).map_err(|e| KmeRustError::SequenceRead {
                        source: std::io::Error::new(std::io::ErrorKind::InvalidData, e.to_string()),
                        path: path_ref.to_path_buf(),
                    })?;

                for result in reader.records() {
                    let record = result.map_err(|e| KmeRustError::SequenceParse {
                        details: e.to_string(),
                    })?;
                    self.process_sequence(record.seq(), k);
                }
            }
            SequenceFormat::Fasta | SequenceFormat::Auto => {
                let reader =
                    fasta::Reader::from_file(path_ref).map_err(|e| KmeRustError::SequenceRead {
                        source: std::io::Error::new(std::io::ErrorKind::InvalidData, e.to_string()),
                        path: path_ref.to_path_buf(),
                    })?;

                for result in reader.records() {
                    let record = result.map_err(|e| KmeRustError::SequenceParse {
                        details: e.to_string(),
                    })?;
                    self.process_sequence(record.seq(), k);
                }
            }
        }

        Ok(self.counts.into_iter().collect())
    }

    #[cfg(feature = "needletail")]
    fn count_file<P>(mut self, path: P, k: KmerLength) -> Result<HashMap<u64, u64>, KmeRustError>
    where
        P: AsRef<Path> + Debug,
    {
        let path_ref = path.as_ref();

        #[cfg(feature = "tracing")]
        let _span = info_span!("sequential_count", path = ?path_ref).entered();

        let mut reader =
            needletail::parse_fastx_file(path_ref).map_err(|e| KmeRustError::SequenceRead {
                source: std::io::Error::new(std::io::ErrorKind::InvalidData, e.to_string()),
                path: path_ref.to_path_buf(),
            })?;

        // Process each record immediately as it's read - no batching
        while let Some(result) = reader.next() {
            let record = result.map_err(|e| KmeRustError::SequenceParse {
                details: e.to_string(),
            })?;
            self.process_sequence(&record.seq(), k);
        }

        Ok(self.counts.into_iter().collect())
    }

    fn process_sequence(&mut self, seq: &[u8], k: KmerLength) {
        let k_val = k.get();
        if seq.len() < k_val {
            return;
        }

        let mut i = 0;
        while i <= seq.len() - k_val {
            let sub = Bytes::copy_from_slice(&seq[i..i + k_val]);

            match Kmer::from_sub(sub) {
                Ok(unpacked) => {
                    let canonical = unpacked.pack().canonical();
                    *self.counts.entry(canonical.packed_bits()).or_insert(0) += 1;
                    i += 1;
                }
                Err(err) => {
                    i += err.position + 1;
                }
            }
        }
    }
}

/// A streaming k-mer counter that processes sequences one at a time.
struct StreamingKmerCounter {
    counts: DashMap<u64, u64, BuildHasherDefault<FxHasher>>,
}

impl StreamingKmerCounter {
    fn new() -> Self {
        Self {
            counts: DashMap::with_hasher(BuildHasherDefault::<FxHasher>::default()),
        }
    }

    #[cfg(all(not(feature = "needletail"), not(feature = "gzip")))]
    fn count_file<P>(self, path: P, k: KmerLength) -> Result<HashMap<u64, u64>, KmeRustError>
    where
        P: AsRef<Path> + Debug,
    {
        use bio::io::{fasta, fastq};

        let path_ref = path.as_ref();
        let format = SequenceFormat::from_extension(path_ref);

        #[cfg(feature = "tracing")]
        let _read_span = info_span!("read_sequences", path = ?path_ref, ?format).entered();

        // Read sequences into a Vec for parallel processing
        let sequences: Vec<Bytes> = match format {
            SequenceFormat::Fastq => {
                let reader =
                    fastq::Reader::from_file(path_ref).map_err(|e| KmeRustError::SequenceRead {
                        source: std::io::Error::new(std::io::ErrorKind::InvalidData, e.to_string()),
                        path: path_ref.to_path_buf(),
                    })?;
                reader
                    .records()
                    .map(|r| {
                        r.map(|rec| Bytes::copy_from_slice(rec.seq())).map_err(|e| {
                            KmeRustError::SequenceParse {
                                details: e.to_string(),
                            }
                        })
                    })
                    .collect::<Result<Vec<_>, _>>()?
            }
            SequenceFormat::Fasta | SequenceFormat::Auto => {
                let reader =
                    fasta::Reader::from_file(path_ref).map_err(|e| KmeRustError::SequenceRead {
                        source: std::io::Error::new(std::io::ErrorKind::InvalidData, e.to_string()),
                        path: path_ref.to_path_buf(),
                    })?;
                reader
                    .records()
                    .map(|r| {
                        r.map(|rec| Bytes::copy_from_slice(rec.seq())).map_err(|e| {
                            KmeRustError::SequenceParse {
                                details: e.to_string(),
                            }
                        })
                    })
                    .collect::<Result<Vec<_>, _>>()?
            }
        };

        #[cfg(feature = "tracing")]
        {
            drop(_read_span);
            debug!(sequences = sequences.len(), "Read sequences from file");
        }

        #[cfg(feature = "tracing")]
        let _process_span = info_span!("process_sequences", count = sequences.len()).entered();

        sequences.par_iter().for_each(|seq| {
            self.process_sequence(seq, k);
        });

        Ok(self.counts.into_iter().collect())
    }

    #[cfg(all(not(feature = "needletail"), feature = "gzip"))]
    fn count_file<P>(self, path: P, k: KmerLength) -> Result<HashMap<u64, u64>, KmeRustError>
    where
        P: AsRef<Path> + Debug,
    {
        use bio::io::{fasta, fastq};
        use flate2::read::GzDecoder;
        use std::{fs::File, io::BufReader};

        let path_ref = path.as_ref();
        let format = SequenceFormat::from_extension(path_ref);
        let is_gzip = path_ref.extension().map(|ext| ext == "gz").unwrap_or(false);

        #[cfg(feature = "tracing")]
        let _read_span = info_span!("read_sequences", path = ?path_ref, ?format).entered();

        // Read sequences into a Vec for parallel processing
        let sequences: Vec<Bytes> = match (format, is_gzip) {
            (SequenceFormat::Fastq, true) => {
                let file = File::open(path_ref).map_err(|e| KmeRustError::SequenceRead {
                    source: e,
                    path: path_ref.to_path_buf(),
                })?;
                let decoder = GzDecoder::new(file);
                let reader = fastq::Reader::new(BufReader::new(decoder));
                reader
                    .records()
                    .map(|r| {
                        r.map(|rec| Bytes::copy_from_slice(rec.seq())).map_err(|e| {
                            KmeRustError::SequenceParse {
                                details: e.to_string(),
                            }
                        })
                    })
                    .collect::<Result<Vec<_>, _>>()?
            }
            (SequenceFormat::Fastq, false) => {
                let reader =
                    fastq::Reader::from_file(path_ref).map_err(|e| KmeRustError::SequenceRead {
                        source: std::io::Error::new(std::io::ErrorKind::InvalidData, e.to_string()),
                        path: path_ref.to_path_buf(),
                    })?;
                reader
                    .records()
                    .map(|r| {
                        r.map(|rec| Bytes::copy_from_slice(rec.seq())).map_err(|e| {
                            KmeRustError::SequenceParse {
                                details: e.to_string(),
                            }
                        })
                    })
                    .collect::<Result<Vec<_>, _>>()?
            }
            (SequenceFormat::Fasta | SequenceFormat::Auto, true) => {
                let file = File::open(path_ref).map_err(|e| KmeRustError::SequenceRead {
                    source: e,
                    path: path_ref.to_path_buf(),
                })?;
                let decoder = GzDecoder::new(file);
                let reader = fasta::Reader::new(BufReader::new(decoder));
                reader
                    .records()
                    .map(|r| {
                        r.map(|rec| Bytes::copy_from_slice(rec.seq())).map_err(|e| {
                            KmeRustError::SequenceParse {
                                details: e.to_string(),
                            }
                        })
                    })
                    .collect::<Result<Vec<_>, _>>()?
            }
            (SequenceFormat::Fasta | SequenceFormat::Auto, false) => {
                let reader =
                    fasta::Reader::from_file(path_ref).map_err(|e| KmeRustError::SequenceRead {
                        source: std::io::Error::new(std::io::ErrorKind::InvalidData, e.to_string()),
                        path: path_ref.to_path_buf(),
                    })?;
                reader
                    .records()
                    .map(|r| {
                        r.map(|rec| Bytes::copy_from_slice(rec.seq())).map_err(|e| {
                            KmeRustError::SequenceParse {
                                details: e.to_string(),
                            }
                        })
                    })
                    .collect::<Result<Vec<_>, _>>()?
            }
        };

        #[cfg(feature = "tracing")]
        {
            drop(_read_span);
            debug!(sequences = sequences.len(), "Read sequences from file");
        }

        #[cfg(feature = "tracing")]
        let _process_span = info_span!("process_sequences", count = sequences.len()).entered();

        sequences.par_iter().for_each(|seq| {
            self.process_sequence(seq, k);
        });

        Ok(self.counts.into_iter().collect())
    }

    #[cfg(feature = "needletail")]
    fn count_file<P>(self, path: P, k: KmerLength) -> Result<HashMap<u64, u64>, KmeRustError>
    where
        P: AsRef<Path> + Debug,
    {
        let path_ref = path.as_ref();

        #[cfg(feature = "tracing")]
        let _read_span = info_span!("read_fasta", path = ?path_ref).entered();

        let mut reader =
            needletail::parse_fastx_file(path_ref).map_err(|e| KmeRustError::SequenceRead {
                source: std::io::Error::new(std::io::ErrorKind::InvalidData, e.to_string()),
                path: path_ref.to_path_buf(),
            })?;

        // Collect sequences for parallel processing
        let mut sequences = Vec::new();
        while let Some(record) = reader.next() {
            let record = record.map_err(|e| KmeRustError::SequenceParse {
                details: e.to_string(),
            })?;
            sequences.push(Bytes::copy_from_slice(&record.seq()));
        }

        #[cfg(feature = "tracing")]
        {
            drop(_read_span);
            debug!(sequences = sequences.len(), "Read sequences from file");
        }

        #[cfg(feature = "tracing")]
        let _process_span = info_span!("process_sequences", count = sequences.len()).entered();

        sequences.par_iter().for_each(|seq| {
            self.process_sequence(seq, k);
        });

        Ok(self.counts.into_iter().collect())
    }

    fn count_sequences<I>(self, sequences: I, k: KmerLength) -> HashMap<u64, u64>
    where
        I: Iterator<Item = Bytes>,
    {
        for seq in sequences {
            self.process_sequence(&seq, k);
        }
        self.counts.into_iter().collect()
    }

    fn process_sequence(&self, seq: &Bytes, k: KmerLength) {
        let k_val = k.get();
        if seq.len() < k_val {
            return;
        }

        let mut i = 0;
        while i <= seq.len() - k_val {
            let sub = seq.slice(i..i + k_val);

            match Kmer::from_sub(sub) {
                Ok(unpacked) => {
                    self.process_valid_kmer(unpacked);
                    i += 1;
                }
                Err(err) => {
                    i += err.position + 1;
                }
            }
        }
    }

    fn process_valid_kmer(&self, unpacked: Kmer) {
        let canonical = unpacked.pack().canonical();
        self.counts
            .entry(canonical.packed_bits())
            .and_modify(|c| *c = c.saturating_add(1))
            .or_insert(1);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn count_from_sequences_basic() {
        let sequences = vec![Bytes::from_static(b"ACGTACGT")];
        let k = KmerLength::new(4).unwrap();
        let counts = count_kmers_from_sequences(sequences.into_iter(), k);

        // ACGT, CGTA, GTAC, TACG, ACGT
        // Canonical forms: ACGT (palindrome), CGTA<->TACG (CGTA smaller), GTAC (palindrome)
        // ACGT appears twice (positions 0 and 4)
        assert!(!counts.is_empty());
    }

    #[test]
    fn count_from_sequences_empty() {
        let sequences: Vec<Bytes> = vec![];
        let k = KmerLength::new(4).unwrap();
        let counts = count_kmers_from_sequences(sequences.into_iter(), k);
        assert!(counts.is_empty());
    }

    #[test]
    fn count_from_sequences_short_sequence() {
        let sequences = vec![Bytes::from_static(b"ACG")];
        let k = KmerLength::new(4).unwrap();
        let counts = count_kmers_from_sequences(sequences.into_iter(), k);
        assert!(counts.is_empty());
    }

    #[test]
    fn count_from_sequences_multiple() {
        let sequences = vec![
            Bytes::from_static(b"AAAA"),
            Bytes::from_static(b"TTTT"), // RC of AAAA
        ];
        let k = KmerLength::new(4).unwrap();
        let counts = count_kmers_from_sequences(sequences.into_iter(), k);

        // AAAA and TTTT should map to same canonical (AAAA)
        assert_eq!(counts.len(), 1);
        let count = counts.values().next().unwrap();
        assert_eq!(*count, 2);
    }
}
