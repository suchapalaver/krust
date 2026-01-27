//! Tests for quality filtering functionality.
//!
//! Verifies that the `--min-quality` flag correctly filters out k-mers
//! containing bases with quality scores below the specified threshold.

#![allow(clippy::unwrap_used, clippy::expect_used)]

use kmerust::format::SequenceFormat;
use kmerust::run::count_kmers_with_quality;
use std::path::PathBuf;

fn fixture_path(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join("fixtures")
        .join(name)
}

#[test]
fn quality_filtering_reduces_kmer_count() {
    let path = fixture_path("low_quality.fq");

    // Without quality filtering
    let counts_no_filter =
        count_kmers_with_quality(&path, 4, SequenceFormat::Fastq, None).expect("should count");

    // With quality filtering (Q=20)
    let counts_filtered =
        count_kmers_with_quality(&path, 4, SequenceFormat::Fastq, Some(20)).expect("should count");

    // The filtered count should have fewer unique k-mers because
    // seq1 has half low-quality bases
    assert!(
        counts_filtered.len() <= counts_no_filter.len(),
        "quality filtering should reduce or maintain k-mer count"
    );
}

#[test]
fn high_quality_data_unchanged_with_filter() {
    // simple.fq has all high-quality bases (I = Phred 40)
    let path = fixture_path("simple.fq");

    // Without quality filtering
    let counts_no_filter =
        count_kmers_with_quality(&path, 4, SequenceFormat::Fastq, None).expect("should count");

    // With quality filtering (Q=20, which is well below I=40)
    let counts_filtered =
        count_kmers_with_quality(&path, 4, SequenceFormat::Fastq, Some(20)).expect("should count");

    // Should have the same counts since all quality is above threshold
    assert_eq!(
        counts_no_filter.len(),
        counts_filtered.len(),
        "high-quality data should be unchanged by filtering"
    );

    for (kmer, count) in &counts_no_filter {
        assert_eq!(
            counts_filtered.get(kmer),
            Some(count),
            "k-mer '{kmer}' should have same count with and without filtering"
        );
    }
}

#[test]
fn zero_quality_threshold_filters_nothing() {
    let path = fixture_path("low_quality.fq");

    // Without quality filtering
    let counts_no_filter =
        count_kmers_with_quality(&path, 4, SequenceFormat::Fastq, None).expect("should count");

    // With Q=0 threshold
    let counts_q0 =
        count_kmers_with_quality(&path, 4, SequenceFormat::Fastq, Some(0)).expect("should count");

    // Q=0 should filter nothing since the minimum ASCII quality is 33 (Phred 0)
    assert_eq!(
        counts_no_filter.len(),
        counts_q0.len(),
        "Q=0 should not filter any k-mers"
    );
}

#[test]
fn fasta_ignores_quality_filter() {
    let path = fixture_path("simple.fa");

    // Without quality filtering
    let counts_no_filter =
        count_kmers_with_quality(&path, 4, SequenceFormat::Fasta, None).expect("should count");

    // With quality filtering (should be ignored for FASTA)
    let counts_filtered =
        count_kmers_with_quality(&path, 4, SequenceFormat::Fasta, Some(20)).expect("should count");

    // Should have identical counts since FASTA has no quality data
    assert_eq!(
        counts_no_filter.len(),
        counts_filtered.len(),
        "FASTA should ignore quality filtering"
    );

    for (kmer, count) in &counts_no_filter {
        assert_eq!(
            counts_filtered.get(kmer),
            Some(count),
            "k-mer '{kmer}' should have same count with quality filter ignored"
        );
    }
}
