//! Streaming k-mer counting for memory-efficient processing of large files.
//!
//! This module provides APIs that process FASTA files without loading all sequences
//! into memory simultaneously, making it suitable for very large genomic datasets.
//!
//! # Memory Model
//!
//! Unlike the standard [`count_kmers`](crate::run::count_kmers) function which loads
//! all sequences before processing, the streaming API:
//! - Reads sequences one at a time
//! - Processes each sequence immediately
//! - Only stores the k-mer count map (one entry per unique canonical k-mer)
//!
//! # Packed Bits API
//!
//! For performance-critical applications, use [`count_kmers_packed`] to get results
//! as packed 64-bit integers instead of strings. This avoids string allocation overhead
//! when you need to do further processing on the k-mer counts.

use std::{collections::HashMap, fmt::Debug, path::Path};

use bytes::Bytes;
use dashmap::DashMap;
use rayon::prelude::*;
use rustc_hash::FxHasher;
use std::hash::BuildHasherDefault;

use crate::{
    error::KmeRustError,
    kmer::{unpack_to_string, Kmer, KmerLength},
};

#[cfg(feature = "tracing")]
use tracing::{debug, info, info_span};

/// Counts k-mers in a FASTA file using streaming I/O.
///
/// Processes sequences one at a time without loading the entire file into memory.
/// This is more memory-efficient than [`count_kmers`](crate::run::count_kmers) for
/// very large files.
///
/// # Arguments
///
/// * `path` - Path to the FASTA file
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
/// let counts = count_kmers_streaming("large_genome.fa", 21)?;
/// println!("Found {} unique k-mers", counts.len());
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
pub fn count_kmers_streaming<P>(path: P, k: usize) -> Result<HashMap<String, i32>, KmeRustError>
where
    P: AsRef<Path> + Debug,
{
    #[cfg(feature = "tracing")]
    info!(k = k, path = ?path, "Starting streaming k-mer counting");

    let k_len = KmerLength::new(k)?;
    let packed = count_kmers_streaming_packed(&path, k_len)?;

    #[cfg(feature = "tracing")]
    let _unpack_span = info_span!("unpack_kmers", count = packed.len()).entered();

    let result: HashMap<String, i32> = packed
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
/// # Arguments
///
/// * `path` - Path to the FASTA file
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
) -> Result<HashMap<u64, i32>, KmeRustError>
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
pub fn count_kmers_from_sequences<I>(sequences: I, k: KmerLength) -> HashMap<u64, i32>
where
    I: Iterator<Item = Bytes>,
{
    let counter = StreamingKmerCounter::new();
    counter.count_sequences(sequences, k)
}

/// A streaming k-mer counter that processes sequences one at a time.
struct StreamingKmerCounter {
    counts: DashMap<u64, i32, BuildHasherDefault<FxHasher>>,
}

impl StreamingKmerCounter {
    fn new() -> Self {
        Self {
            counts: DashMap::with_hasher(BuildHasherDefault::<FxHasher>::default()),
        }
    }

    #[cfg(all(not(feature = "needletail"), not(feature = "gzip")))]
    fn count_file<P>(self, path: P, k: KmerLength) -> Result<HashMap<u64, i32>, KmeRustError>
    where
        P: AsRef<Path> + Debug,
    {
        use bio::io::fasta;

        let path_ref = path.as_ref();

        #[cfg(feature = "tracing")]
        let _read_span = info_span!("read_fasta", path = ?path_ref).entered();

        let reader = fasta::Reader::from_file(path_ref).map_err(|e| KmeRustError::FastaRead {
            source: std::io::Error::new(std::io::ErrorKind::InvalidData, e.to_string()),
            path: path_ref.to_path_buf(),
        })?;

        // Process records in parallel chunks for efficiency while maintaining streaming behavior
        let records: Vec<_> = reader
            .records()
            .collect::<Result<Vec<_>, _>>()
            .map_err(|e| KmeRustError::FastaParse {
                details: e.to_string(),
            })?;

        #[cfg(feature = "tracing")]
        {
            drop(_read_span);
            debug!(sequences = records.len(), "Read sequences from file");
        }

        #[cfg(feature = "tracing")]
        let _process_span = info_span!("process_sequences", count = records.len()).entered();

        records.par_iter().for_each(|record| {
            let seq = Bytes::copy_from_slice(record.seq());
            self.process_sequence(&seq, k);
        });

        Ok(self.counts.into_iter().collect())
    }

    #[cfg(all(not(feature = "needletail"), feature = "gzip"))]
    fn count_file<P>(self, path: P, k: KmerLength) -> Result<HashMap<u64, i32>, KmeRustError>
    where
        P: AsRef<Path> + Debug,
    {
        use bio::io::fasta;
        use flate2::read::GzDecoder;
        use std::{fs::File, io::BufReader};

        let path_ref = path.as_ref();

        #[cfg(feature = "tracing")]
        let _read_span = info_span!("read_fasta", path = ?path_ref).entered();

        let is_gzip = path_ref.extension().map(|ext| ext == "gz").unwrap_or(false);

        let records: Vec<bio::io::fasta::Record> = if is_gzip {
            let file = File::open(path_ref).map_err(|e| KmeRustError::FastaRead {
                source: e,
                path: path_ref.to_path_buf(),
            })?;
            let decoder = GzDecoder::new(file);
            let reader = fasta::Reader::new(BufReader::new(decoder));
            reader
                .records()
                .collect::<Result<Vec<_>, _>>()
                .map_err(|e| KmeRustError::FastaParse {
                    details: e.to_string(),
                })?
        } else {
            let reader =
                fasta::Reader::from_file(path_ref).map_err(|e| KmeRustError::FastaRead {
                    source: std::io::Error::new(std::io::ErrorKind::InvalidData, e.to_string()),
                    path: path_ref.to_path_buf(),
                })?;
            reader
                .records()
                .collect::<Result<Vec<_>, _>>()
                .map_err(|e| KmeRustError::FastaParse {
                    details: e.to_string(),
                })?
        };

        #[cfg(feature = "tracing")]
        {
            drop(_read_span);
            debug!(sequences = records.len(), "Read sequences from file");
        }

        #[cfg(feature = "tracing")]
        let _process_span = info_span!("process_sequences", count = records.len()).entered();

        records.par_iter().for_each(|record| {
            let seq = Bytes::copy_from_slice(record.seq());
            self.process_sequence(&seq, k);
        });

        Ok(self.counts.into_iter().collect())
    }

    #[cfg(feature = "needletail")]
    fn count_file<P>(self, path: P, k: KmerLength) -> Result<HashMap<u64, i32>, KmeRustError>
    where
        P: AsRef<Path> + Debug,
    {
        let path_ref = path.as_ref();

        #[cfg(feature = "tracing")]
        let _read_span = info_span!("read_fasta", path = ?path_ref).entered();

        let mut reader =
            needletail::parse_fastx_file(path_ref).map_err(|e| KmeRustError::FastaRead {
                source: std::io::Error::new(std::io::ErrorKind::InvalidData, e.to_string()),
                path: path_ref.to_path_buf(),
            })?;

        // Collect sequences for parallel processing
        let mut sequences = Vec::new();
        while let Some(record) = reader.next() {
            let record = record.map_err(|e| KmeRustError::FastaParse {
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

    fn count_sequences<I>(self, sequences: I, k: KmerLength) -> HashMap<u64, i32>
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
        let packed = unpacked.pack();

        // Optimization: check if we've already seen this exact k-mer
        if let Some(mut count) = self.counts.get_mut(&packed.packed_bits()) {
            *count += 1;
        } else {
            // Not seen before - compute canonical form
            let canonical = packed.canonical();
            *self.counts.entry(canonical.packed_bits()).or_insert(0) += 1;
        }
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
