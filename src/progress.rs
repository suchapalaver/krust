//! Progress tracking for k-mer counting operations.
//!
//! This module provides thread-safe progress reporting via callbacks,
//! allowing callers to monitor the progress of long-running k-mer counting operations.
//!
//! # Example
//!
//! ```rust,no_run
//! use kmerust::progress::Progress;
//! use kmerust::run::count_kmers_with_progress;
//!
//! let counts = count_kmers_with_progress("genome.fa", 21, |progress| {
//!     println!(
//!         "Processed {} sequences ({} bases)",
//!         progress.sequences_processed,
//!         progress.bases_processed
//!     );
//! })?;
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```

use std::sync::atomic::{AtomicU64, Ordering};

/// Progress snapshot during k-mer counting.
#[derive(Debug, Clone, Default)]
pub struct Progress {
    /// Number of sequences processed so far.
    pub sequences_processed: u64,
    /// Total number of bases processed so far.
    pub bases_processed: u64,
}

/// Thread-safe progress tracker using atomic counters.
///
/// This struct maintains atomic counters that can be safely updated from
/// multiple threads during parallel k-mer counting.
#[derive(Debug, Default)]
pub struct ProgressTracker {
    sequences: AtomicU64,
    bases: AtomicU64,
}

impl ProgressTracker {
    /// Create a new progress tracker with zero counts.
    #[must_use]
    pub const fn new() -> Self {
        Self {
            sequences: AtomicU64::new(0),
            bases: AtomicU64::new(0),
        }
    }

    /// Record that a sequence has been processed.
    ///
    /// This method is thread-safe and can be called from multiple threads.
    ///
    /// # Arguments
    ///
    /// * `bases` - The number of bases in the processed sequence.
    pub fn record_sequence(&self, bases: u64) {
        self.sequences.fetch_add(1, Ordering::Relaxed);
        self.bases.fetch_add(bases, Ordering::Relaxed);
    }

    /// Get a snapshot of the current progress.
    ///
    /// The returned values represent the state at a point in time and may
    /// change immediately after this call returns.
    pub fn snapshot(&self) -> Progress {
        Progress {
            sequences_processed: self.sequences.load(Ordering::Relaxed),
            bases_processed: self.bases.load(Ordering::Relaxed),
        }
    }

    /// Reset all counters to zero.
    pub fn reset(&self) {
        self.sequences.store(0, Ordering::Relaxed);
        self.bases.store(0, Ordering::Relaxed);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tracker_starts_at_zero() {
        let tracker = ProgressTracker::new();
        let progress = tracker.snapshot();
        assert_eq!(progress.sequences_processed, 0);
        assert_eq!(progress.bases_processed, 0);
    }

    #[test]
    fn tracker_records_sequence() {
        let tracker = ProgressTracker::new();
        tracker.record_sequence(100);
        tracker.record_sequence(50);

        let progress = tracker.snapshot();
        assert_eq!(progress.sequences_processed, 2);
        assert_eq!(progress.bases_processed, 150);
    }

    #[test]
    fn tracker_reset() {
        let tracker = ProgressTracker::new();
        tracker.record_sequence(100);
        tracker.reset();

        let progress = tracker.snapshot();
        assert_eq!(progress.sequences_processed, 0);
        assert_eq!(progress.bases_processed, 0);
    }
}
