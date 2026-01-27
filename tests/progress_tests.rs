//! Tests for progress reporting during k-mer counting.

#![allow(clippy::unwrap_used, clippy::expect_used)]

use kmerust::progress::{Progress, ProgressTracker};
use kmerust::run::count_kmers_with_progress;
use std::path::PathBuf;
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::Arc;

fn fixture_path(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join("fixtures")
        .join(name)
}

#[test]
fn progress_tracker_records_sequences() {
    let tracker = ProgressTracker::new();

    tracker.record_sequence(100);
    tracker.record_sequence(200);
    tracker.record_sequence(50);

    let progress = tracker.snapshot();
    assert_eq!(progress.sequences_processed, 3);
    assert_eq!(progress.bases_processed, 350);
}

#[test]
fn progress_tracker_reset() {
    let tracker = ProgressTracker::new();
    tracker.record_sequence(100);

    tracker.reset();

    let progress = tracker.snapshot();
    assert_eq!(progress.sequences_processed, 0);
    assert_eq!(progress.bases_processed, 0);
}

#[test]
fn count_kmers_with_progress_invokes_callback() {
    let path = fixture_path("simple.fa");
    let callback_count = Arc::new(AtomicU64::new(0));
    let callback_count_clone = Arc::clone(&callback_count);

    let counts = count_kmers_with_progress(&path, 4, move |_progress| {
        callback_count_clone.fetch_add(1, Ordering::SeqCst);
    })
    .expect("should count k-mers");

    // Should have invoked callback at least once
    assert!(
        callback_count.load(Ordering::SeqCst) > 0,
        "progress callback should be invoked"
    );

    // Should still produce correct results
    assert!(!counts.is_empty(), "should find k-mers");
}

#[test]
fn count_kmers_with_progress_tracks_sequences() {
    let path = fixture_path("simple.fa");
    let max_sequences = Arc::new(AtomicU64::new(0));
    let max_sequences_clone = Arc::clone(&max_sequences);

    let _counts = count_kmers_with_progress(&path, 4, move |progress: Progress| {
        // Track the maximum number of sequences seen
        max_sequences_clone.fetch_max(progress.sequences_processed, Ordering::SeqCst);
    })
    .expect("should count k-mers");

    // simple.fa has 2 sequences
    assert_eq!(
        max_sequences.load(Ordering::SeqCst),
        2,
        "should have processed 2 sequences"
    );
}

#[test]
fn count_kmers_with_progress_produces_same_results_as_regular() {
    use kmerust::run::count_kmers;

    let path = fixture_path("simple.fa");

    let regular_counts = count_kmers(&path, 4).expect("regular count should work");
    let progress_counts =
        count_kmers_with_progress(&path, 4, |_| {}).expect("progress count should work");

    assert_eq!(
        regular_counts.len(),
        progress_counts.len(),
        "should have same number of k-mers"
    );

    for (kmer, count) in &regular_counts {
        assert_eq!(
            progress_counts.get(kmer),
            Some(count),
            "k-mer '{kmer}' should have same count"
        );
    }
}
