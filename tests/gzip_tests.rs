//! Tests for gzip compressed input support.

#![cfg(feature = "gzip")]

use kmerust::run::count_kmers;
use kmerust::streaming::count_kmers_streaming;
use std::path::PathBuf;

fn fixture_path(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join("fixtures")
        .join(name)
}

#[test]
fn count_kmers_from_gzip_file() {
    let path = fixture_path("simple.fa.gz");
    let counts = count_kmers(&path, 4).expect("should count k-mers from gzipped file");

    // Should have some k-mers
    assert!(!counts.is_empty(), "should find k-mers in gzipped file");
}

#[test]
fn gzip_and_plain_produce_same_results() {
    let plain_path = fixture_path("simple.fa");
    let gzip_path = fixture_path("simple.fa.gz");

    let plain_counts = count_kmers(&plain_path, 4).expect("should count plain file");
    let gzip_counts = count_kmers(&gzip_path, 4).expect("should count gzipped file");

    assert_eq!(
        plain_counts.len(),
        gzip_counts.len(),
        "should have same number of unique k-mers"
    );

    for (kmer, count) in &plain_counts {
        assert_eq!(
            gzip_counts.get(kmer),
            Some(count),
            "k-mer '{kmer}' should have same count in both files"
        );
    }
}

#[test]
fn streaming_count_from_gzip_file() {
    let path = fixture_path("simple.fa.gz");
    let counts = count_kmers_streaming(&path, 4).expect("should stream k-mers from gzipped file");

    // Should have some k-mers
    assert!(!counts.is_empty(), "should find k-mers in gzipped file");
}

#[test]
fn streaming_gzip_and_plain_produce_same_results() {
    let plain_path = fixture_path("simple.fa");
    let gzip_path = fixture_path("simple.fa.gz");

    let plain_counts = count_kmers_streaming(&plain_path, 4).expect("should stream plain file");
    let gzip_counts = count_kmers_streaming(&gzip_path, 4).expect("should stream gzipped file");

    assert_eq!(
        plain_counts.len(),
        gzip_counts.len(),
        "should have same number of unique k-mers"
    );

    for (kmer, count) in &plain_counts {
        assert_eq!(
            gzip_counts.get(kmer),
            Some(count),
            "k-mer '{kmer}' should have same count in both files"
        );
    }
}
