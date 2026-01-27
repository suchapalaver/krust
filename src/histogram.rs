//! K-mer frequency histogram computation.
//!
//! This module provides types and functions for computing k-mer frequency histograms
//! (count of counts). K-mer histograms are fundamental for genome size estimation,
//! error detection, and heterozygosity analysis.
//!
//! # Example
//!
//! ```rust
//! use std::collections::HashMap;
//! use kmerust::histogram::{compute_histogram, KmerHistogram};
//!
//! // Suppose we have k-mer counts
//! let counts: HashMap<String, u64> = [
//!     ("ACG".to_string(), 1),
//!     ("CGT".to_string(), 1),
//!     ("GTA".to_string(), 2),
//!     ("TAC".to_string(), 2),
//! ].into();
//!
//! let histogram = compute_histogram(&counts);
//!
//! // 2 k-mers appear once, 2 k-mers appear twice
//! assert_eq!(histogram.get(&1), Some(&2));
//! assert_eq!(histogram.get(&2), Some(&2));
//! ```

use std::collections::{BTreeMap, HashMap};

/// K-mer frequency histogram: maps count -> number of distinct k-mers with that count.
///
/// Uses `BTreeMap` for sorted iteration (counts in ascending order).
pub type KmerHistogram = BTreeMap<u64, u64>;

/// Summary statistics for a k-mer histogram.
///
/// These statistics are useful for genome analysis:
/// - `total_kmers`: Total k-mer occurrences (sum of all counts)
/// - `distinct_kmers`: Number of unique k-mers
/// - `mode_count`: The count value that appears most frequently
/// - `mode_frequency`: How many k-mers have the mode count
/// - `mean_count`: Average count per unique k-mer
#[derive(Debug, Clone, PartialEq)]
pub struct HistogramStats {
    /// Total k-mer occurrences (sum of all k-mer counts).
    pub total_kmers: u64,
    /// Number of unique k-mers.
    pub distinct_kmers: u64,
    /// The count value that appears most frequently (mode of the distribution).
    pub mode_count: u64,
    /// Number of k-mers that have the mode count.
    pub mode_frequency: u64,
    /// Average k-mer count (`total_kmers` / `distinct_kmers`).
    pub mean_count: f64,
}

/// Computes a histogram from k-mer counts.
///
/// Given a mapping of k-mer strings to their counts, returns a histogram
/// showing how many k-mers have each count value.
///
/// # Arguments
///
/// * `counts` - A `HashMap` mapping k-mer strings to their counts
///
/// # Returns
///
/// A `KmerHistogram` (`BTreeMap`) mapping count -> frequency.
///
/// # Example
///
/// ```rust
/// use std::collections::HashMap;
/// use kmerust::histogram::compute_histogram;
///
/// let counts: HashMap<String, u64> = [
///     ("ACG".to_string(), 5),
///     ("CGT".to_string(), 5),
///     ("GTA".to_string(), 10),
/// ].into();
///
/// let hist = compute_histogram(&counts);
/// assert_eq!(hist.get(&5), Some(&2));  // 2 k-mers with count 5
/// assert_eq!(hist.get(&10), Some(&1)); // 1 k-mer with count 10
/// ```
#[must_use]
#[allow(clippy::implicit_hasher)]
pub fn compute_histogram(counts: &HashMap<String, u64>) -> KmerHistogram {
    let mut histogram = BTreeMap::new();
    for &count in counts.values() {
        *histogram.entry(count).or_insert(0) += 1;
    }
    histogram
}

/// Computes a histogram from packed k-mer counts.
///
/// This is more memory efficient when working with packed (bit-encoded) k-mer
/// representations directly.
///
/// # Arguments
///
/// * `counts` - A `HashMap` mapping packed k-mer bits to their counts
///
/// # Returns
///
/// A `KmerHistogram` (`BTreeMap`) mapping count -> frequency.
#[must_use]
#[allow(clippy::implicit_hasher)]
pub fn compute_histogram_packed(counts: &HashMap<u64, u64>) -> KmerHistogram {
    let mut histogram = BTreeMap::new();
    for &count in counts.values() {
        *histogram.entry(count).or_insert(0) += 1;
    }
    histogram
}

/// Computes summary statistics for a k-mer histogram.
///
/// # Arguments
///
/// * `histogram` - A reference to a `KmerHistogram`
///
/// # Returns
///
/// A `HistogramStats` struct with summary statistics.
///
/// # Example
///
/// ```rust
/// use std::collections::HashMap;
/// use kmerust::histogram::{compute_histogram, histogram_stats};
///
/// let counts: HashMap<String, u64> = [
///     ("ACG".to_string(), 1),
///     ("CGT".to_string(), 1),
///     ("GTA".to_string(), 2),
///     ("TAC".to_string(), 2),
/// ].into();
///
/// let hist = compute_histogram(&counts);
/// let stats = histogram_stats(&hist);
///
/// assert_eq!(stats.distinct_kmers, 4);
/// assert_eq!(stats.total_kmers, 6); // 1+1+2+2
/// ```
#[must_use]
pub fn histogram_stats(histogram: &KmerHistogram) -> HistogramStats {
    let distinct: u64 = histogram.values().sum();
    let total: u64 = histogram.iter().map(|(c, f)| c * f).sum();

    let (mode_count, mode_frequency) = histogram
        .iter()
        .max_by_key(|(_, f)| *f)
        .map_or((0, 0), |(&c, &f)| (c, f));

    HistogramStats {
        total_kmers: total,
        distinct_kmers: distinct,
        mode_count,
        mode_frequency,
        #[allow(clippy::cast_precision_loss)]
        mean_count: if distinct > 0 {
            total as f64 / distinct as f64
        } else {
            0.0
        },
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn histogram_basic() {
        let counts: HashMap<String, u64> = [
            ("ACG".to_string(), 1),
            ("CGT".to_string(), 1),
            ("GTA".to_string(), 2),
            ("TAC".to_string(), 2),
        ]
        .into();

        let hist = compute_histogram(&counts);

        // 2 k-mers appear once
        assert_eq!(hist.get(&1), Some(&2));
        // 2 k-mers appear twice
        assert_eq!(hist.get(&2), Some(&2));
        // No k-mers appear 3 times
        assert_eq!(hist.get(&3), None);
    }

    #[test]
    fn histogram_single_kmer() {
        let counts: HashMap<String, u64> = [("ACGT".to_string(), 100)].into();

        let hist = compute_histogram(&counts);

        assert_eq!(hist.len(), 1);
        assert_eq!(hist.get(&100), Some(&1));
    }

    #[test]
    fn histogram_empty() {
        let counts: HashMap<String, u64> = HashMap::new();
        let hist = compute_histogram(&counts);
        assert!(hist.is_empty());
    }

    #[test]
    fn histogram_packed() {
        let counts: HashMap<u64, u64> = [
            (0b0001, 5),  // Some packed k-mer with count 5
            (0b0010, 5),  // Another packed k-mer with count 5
            (0b0011, 10), // Another with count 10
        ]
        .into();

        let hist = compute_histogram_packed(&counts);

        assert_eq!(hist.get(&5), Some(&2));
        assert_eq!(hist.get(&10), Some(&1));
    }

    #[test]
    fn histogram_stats_basic() {
        let counts: HashMap<String, u64> = [
            ("ACG".to_string(), 1),
            ("CGT".to_string(), 1),
            ("GTA".to_string(), 2),
            ("TAC".to_string(), 2),
        ]
        .into();

        let hist = compute_histogram(&counts);
        let stats = histogram_stats(&hist);

        assert_eq!(stats.distinct_kmers, 4);
        assert_eq!(stats.total_kmers, 6); // 1+1+2+2
                                          // Both count 1 and count 2 have frequency 2, mode is one of them
        assert!(stats.mode_frequency == 2);
        assert!((stats.mean_count - 1.5).abs() < f64::EPSILON);
    }

    #[test]
    fn histogram_stats_empty() {
        let hist = KmerHistogram::new();
        let stats = histogram_stats(&hist);

        assert_eq!(stats.distinct_kmers, 0);
        assert_eq!(stats.total_kmers, 0);
        assert_eq!(stats.mode_count, 0);
        assert_eq!(stats.mode_frequency, 0);
        assert!((stats.mean_count - 0.0).abs() < f64::EPSILON);
    }

    #[test]
    fn histogram_stats_single_kmer() {
        let counts: HashMap<String, u64> = [("ACGT".to_string(), 42)].into();

        let hist = compute_histogram(&counts);
        let stats = histogram_stats(&hist);

        assert_eq!(stats.distinct_kmers, 1);
        assert_eq!(stats.total_kmers, 42);
        assert_eq!(stats.mode_count, 42);
        assert_eq!(stats.mode_frequency, 1);
        assert!((stats.mean_count - 42.0).abs() < f64::EPSILON);
    }

    #[test]
    fn histogram_sorted_keys() {
        let counts: HashMap<String, u64> = [
            ("A".to_string(), 100),
            ("B".to_string(), 1),
            ("C".to_string(), 50),
        ]
        .into();

        let hist = compute_histogram(&counts);
        let keys: Vec<_> = hist.keys().collect();

        // BTreeMap should have sorted keys
        assert_eq!(keys, vec![&1, &50, &100]);
    }
}
