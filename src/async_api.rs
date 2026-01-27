//! Async k-mer counting API using Tokio.
//!
//! This module provides async versions of the k-mer counting functions,
//! suitable for use in async runtimes like Tokio.
//!
//! # Feature Flag
//!
//! This module requires the `async` feature to be enabled:
//!
//! ```toml
//! [dependencies]
//! kmerust = { version = "0.1", features = ["async"] }
//! ```
//!
//! # Example
//!
//! ```rust,no_run
//! use kmerust::async_api::count_kmers_async;
//!
//! #[tokio::main]
//! async fn main() -> Result<(), Box<dyn std::error::Error + Send + Sync>> {
//!     let counts = count_kmers_async("genome.fa", 21).await?;
//!     println!("Found {} unique k-mers", counts.len());
//!     Ok(())
//! }
//! ```

use std::{collections::HashMap, fmt::Debug, path::Path};

use tokio::task;

use crate::{
    error::KmeRustError,
    kmer::{unpack_to_string, KmerLength},
    streaming::count_kmers_streaming_packed,
};

/// Async version of [`count_kmers`](crate::run::count_kmers).
///
/// Spawns the CPU-intensive k-mer counting work on Tokio's blocking thread pool,
/// allowing other async tasks to make progress while counting is in progress.
///
/// # Arguments
///
/// * `path` - Path to the FASTA file
/// * `k` - K-mer length (must be 1-32)
///
/// # Returns
///
/// A `HashMap` mapping k-mer strings to their counts.
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
/// use kmerust::async_api::count_kmers_async;
///
/// #[tokio::main]
/// async fn main() -> Result<(), Box<dyn std::error::Error + Send + Sync>> {
///     let counts = count_kmers_async("genome.fa", 21).await?;
///
///     for (kmer, count) in counts.iter().take(10) {
///         println!("{kmer}: {count}");
///     }
///     Ok(())
/// }
/// ```
pub async fn count_kmers_async<P>(
    path: P,
    k: usize,
) -> Result<HashMap<String, u64>, Box<dyn std::error::Error + Send + Sync>>
where
    P: AsRef<Path> + Debug + Send + 'static,
{
    let k_len = KmerLength::new(k)?;

    // Spawn blocking work on Tokio's thread pool
    let packed = task::spawn_blocking(move || count_kmers_streaming_packed(path, k_len)).await??;

    // Convert packed bits to strings
    let result: HashMap<String, u64> = packed
        .into_iter()
        .map(|(bits, count)| (unpack_to_string(bits, k_len), count))
        .collect();

    Ok(result)
}

/// Async version of [`crate::streaming::count_kmers_streaming_packed`].
///
/// Returns packed 64-bit representations for maximum efficiency.
///
/// # Arguments
///
/// * `path` - Path to the FASTA file
/// * `k` - Validated k-mer length
///
/// # Returns
///
/// A `HashMap` mapping packed k-mer bits to their counts.
///
/// # Example
///
/// ```rust,no_run
/// use kmerust::async_api::count_kmers_packed_async;
/// use kmerust::kmer::KmerLength;
///
/// #[tokio::main]
/// async fn main() -> Result<(), Box<dyn std::error::Error + Send + Sync>> {
///     let k = KmerLength::new(21)?;
///     let counts = count_kmers_packed_async("genome.fa", k).await?;
///
///     for (packed_bits, count) in counts.iter().take(10) {
///         println!("{packed_bits:#x}: {count}");
///     }
///     Ok(())
/// }
/// ```
pub async fn count_kmers_packed_async<P>(
    path: P,
    k: KmerLength,
) -> Result<HashMap<u64, u64>, Box<dyn std::error::Error + Send + Sync>>
where
    P: AsRef<Path> + Debug + Send + 'static,
{
    let result = task::spawn_blocking(move || count_kmers_streaming_packed(path, k)).await??;
    Ok(result)
}

/// Async builder for k-mer counting operations.
///
/// Provides the same fluent API as [`KmerCounter`](crate::builder::KmerCounter)
/// but with async execution.
///
/// # Example
///
/// ```rust,no_run
/// use kmerust::async_api::AsyncKmerCounter;
///
/// #[tokio::main]
/// async fn main() -> Result<(), Box<dyn std::error::Error + Send + Sync>> {
///     let counts = AsyncKmerCounter::new()
///         .k(21)?
///         .min_count(5)
///         .count("genome.fa")
///         .await?;
///
///     println!("Found {} unique k-mers", counts.len());
///     Ok(())
/// }
/// ```
#[derive(Debug, Clone)]
pub struct AsyncKmerCounter {
    k: Option<KmerLength>,
    min_count: u64,
}

impl Default for AsyncKmerCounter {
    fn default() -> Self {
        Self::new()
    }
}

impl AsyncKmerCounter {
    /// Creates a new async k-mer counter builder.
    #[must_use]
    pub const fn new() -> Self {
        Self {
            k: None,
            min_count: 1,
        }
    }

    /// Sets the k-mer length.
    ///
    /// # Errors
    ///
    /// Returns [`KmerLengthError`](crate::error::KmerLengthError) if `k` is outside 1-32.
    pub fn k(mut self, k: usize) -> Result<Self, crate::error::KmerLengthError> {
        self.k = Some(KmerLength::new(k)?);
        Ok(self)
    }

    /// Sets the k-mer length from a pre-validated `KmerLength`.
    #[must_use]
    pub const fn k_validated(mut self, k: KmerLength) -> Self {
        self.k = Some(k);
        self
    }

    /// Sets the minimum count threshold.
    #[must_use]
    pub const fn min_count(mut self, min_count: u64) -> Self {
        self.min_count = min_count;
        self
    }

    /// Counts k-mers asynchronously.
    ///
    /// # Errors
    ///
    /// Returns an error if:
    /// - `k` has not been set
    /// - The file cannot be read or parsed
    pub async fn count<P>(
        &self,
        path: P,
    ) -> Result<HashMap<String, u64>, Box<dyn std::error::Error + Send + Sync>>
    where
        P: AsRef<Path> + Debug + Send + 'static,
    {
        let k = self.k.ok_or(KmeRustError::InvalidKmerLength {
            k: 0,
            min: 1,
            max: 32,
        })?;
        let min_count = self.min_count;

        let counts = count_kmers_async(path, k.get()).await?;

        if min_count > 1 {
            Ok(counts
                .into_iter()
                .filter(|(_, count)| *count >= min_count)
                .collect())
        } else {
            Ok(counts)
        }
    }

    /// Counts k-mers and returns packed representations asynchronously.
    ///
    /// # Errors
    ///
    /// Returns an error if:
    /// - `k` has not been set
    /// - The file cannot be read or parsed
    pub async fn count_packed<P>(
        &self,
        path: P,
    ) -> Result<HashMap<u64, u64>, Box<dyn std::error::Error + Send + Sync>>
    where
        P: AsRef<Path> + Debug + Send + 'static,
    {
        let k = self.k.ok_or(KmeRustError::InvalidKmerLength {
            k: 0,
            min: 1,
            max: 32,
        })?;
        let min_count = self.min_count;

        let counts = count_kmers_packed_async(path, k).await?;

        if min_count > 1 {
            Ok(counts
                .into_iter()
                .filter(|(_, count)| *count >= min_count)
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
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;

    #[test]
    fn async_counter_default() {
        let counter = AsyncKmerCounter::new();
        assert!(counter.get_k().is_none());
        assert_eq!(counter.get_min_count(), 1);
    }

    #[test]
    fn async_counter_k_valid() {
        let counter = AsyncKmerCounter::new().k(21).unwrap();
        assert_eq!(counter.get_k().unwrap().get(), 21);
    }

    #[test]
    fn async_counter_k_invalid() {
        let result = AsyncKmerCounter::new().k(0);
        assert!(result.is_err());

        let result = AsyncKmerCounter::new().k(33);
        assert!(result.is_err());
    }

    #[test]
    fn async_counter_chained() {
        let counter = AsyncKmerCounter::new().k(21).unwrap().min_count(5);

        assert_eq!(counter.get_k().unwrap().get(), 21);
        assert_eq!(counter.get_min_count(), 5);
    }
}
