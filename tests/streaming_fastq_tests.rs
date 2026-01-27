//! Tests for streaming FASTQ support.
//!
//! Verifies that the streaming APIs (`count_kmers_streaming`, `count_kmers_sequential`)
//! correctly handle FASTQ files with automatic format detection.

use kmerust::streaming::{count_kmers_sequential, count_kmers_streaming};
use std::path::PathBuf;

fn fixture_path(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join("fixtures")
        .join(name)
}

#[test]
fn streaming_count_from_fastq() {
    let path = fixture_path("simple.fq");
    let counts = count_kmers_streaming(&path, 4).expect("should stream k-mers from FASTQ file");

    assert!(!counts.is_empty(), "should find k-mers in FASTQ file");
}

#[test]
fn sequential_count_from_fastq() {
    let path = fixture_path("simple.fq");
    let counts =
        count_kmers_sequential(&path, 4).expect("should count k-mers sequentially from FASTQ");

    assert!(!counts.is_empty(), "should find k-mers in FASTQ file");
}

#[test]
fn streaming_fastq_handles_n_bases() {
    let path = fixture_path("with_n.fq");
    let counts = count_kmers_streaming(&path, 4).expect("should handle FASTQ with N bases");

    // Should have k-mers, but N bases cause breaks in k-mer extraction
    assert!(!counts.is_empty(), "should find k-mers despite N bases");
}

#[test]
fn sequential_fastq_handles_n_bases() {
    let path = fixture_path("with_n.fq");
    let counts =
        count_kmers_sequential(&path, 4).expect("should handle FASTQ with N bases sequentially");

    assert!(!counts.is_empty(), "should find k-mers despite N bases");
}

#[test]
fn streaming_fastq_matches_fasta_counts() {
    // Compare FASTQ counts with equivalent FASTA
    let fasta_path = fixture_path("simple.fa");
    let fastq_path = fixture_path("simple.fq");

    let fasta_counts = count_kmers_streaming(&fasta_path, 4).expect("should stream FASTA");
    let fastq_counts = count_kmers_streaming(&fastq_path, 4).expect("should stream FASTQ");

    // The fixtures should have the same sequences
    assert_eq!(
        fasta_counts.len(),
        fastq_counts.len(),
        "FASTA and FASTQ should have same unique k-mer count"
    );

    for (kmer, count) in &fasta_counts {
        assert_eq!(
            fastq_counts.get(kmer),
            Some(count),
            "k-mer '{kmer}' should have same count in both formats"
        );
    }
}

#[test]
fn sequential_fastq_matches_fasta_counts() {
    let fasta_path = fixture_path("simple.fa");
    let fastq_path = fixture_path("simple.fq");

    let fasta_counts =
        count_kmers_sequential(&fasta_path, 4).expect("should count FASTA sequentially");
    let fastq_counts =
        count_kmers_sequential(&fastq_path, 4).expect("should count FASTQ sequentially");

    assert_eq!(
        fasta_counts.len(),
        fastq_counts.len(),
        "FASTA and FASTQ should have same unique k-mer count"
    );

    for (packed_bits, count) in &fasta_counts {
        assert_eq!(
            fastq_counts.get(packed_bits),
            Some(count),
            "packed k-mer should have same count in both formats"
        );
    }
}

#[cfg(feature = "gzip")]
mod gzip_fastq {
    use super::*;

    #[test]
    fn streaming_count_from_gzip_fastq() {
        let path = fixture_path("simple.fq.gz");
        let counts =
            count_kmers_streaming(&path, 4).expect("should stream k-mers from gzipped FASTQ");

        assert!(!counts.is_empty(), "should find k-mers in gzipped FASTQ");
    }

    #[test]
    fn sequential_count_from_gzip_fastq() {
        let path = fixture_path("simple.fq.gz");
        let counts = count_kmers_sequential(&path, 4)
            .expect("should count k-mers sequentially from gzipped FASTQ");

        assert!(!counts.is_empty(), "should find k-mers in gzipped FASTQ");
    }

    #[test]
    fn streaming_gzip_fastq_matches_plain_fastq() {
        let plain_path = fixture_path("simple.fq");
        let gzip_path = fixture_path("simple.fq.gz");

        let plain_counts =
            count_kmers_streaming(&plain_path, 4).expect("should stream plain FASTQ");
        let gzip_counts =
            count_kmers_streaming(&gzip_path, 4).expect("should stream gzipped FASTQ");

        assert_eq!(
            plain_counts.len(),
            gzip_counts.len(),
            "plain and gzipped FASTQ should have same unique k-mer count"
        );

        for (kmer, count) in &plain_counts {
            assert_eq!(
                gzip_counts.get(kmer),
                Some(count),
                "k-mer '{kmer}' should have same count in both plain and gzipped FASTQ"
            );
        }
    }

    #[test]
    fn sequential_gzip_fastq_matches_plain_fastq() {
        let plain_path = fixture_path("simple.fq");
        let gzip_path = fixture_path("simple.fq.gz");

        let plain_counts =
            count_kmers_sequential(&plain_path, 4).expect("should count plain FASTQ");
        let gzip_counts =
            count_kmers_sequential(&gzip_path, 4).expect("should count gzipped FASTQ");

        assert_eq!(
            plain_counts.len(),
            gzip_counts.len(),
            "plain and gzipped FASTQ should have same unique k-mer count"
        );

        for (packed_bits, count) in &plain_counts {
            assert_eq!(
                gzip_counts.get(packed_bits),
                Some(count),
                "packed k-mer should have same count in both plain and gzipped FASTQ"
            );
        }
    }
}
