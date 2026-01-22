//! Tests for memory-mapped I/O support.

#![cfg(feature = "mmap")]

use kmerust::mmap::MmapFasta;
use kmerust::run::{count_kmers, count_kmers_mmap};
use std::io::Write;
use std::path::PathBuf;
use tempfile::NamedTempFile;

fn fixture_path(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join("fixtures")
        .join(name)
}

#[test]
fn mmap_fasta_open_and_read() {
    let path = fixture_path("simple.fa");
    let mmap = MmapFasta::open(&path).expect("should open file");

    assert!(!mmap.is_empty(), "file should not be empty");
    assert!(
        mmap.as_bytes().starts_with(b">seq1"),
        "should start with fasta header"
    );
}

#[test]
fn mmap_fasta_len() {
    let mut temp = NamedTempFile::new().unwrap();
    write!(temp, "ACGT").unwrap();
    temp.flush().unwrap();

    let mmap = MmapFasta::open(temp.path()).expect("should open temp file");
    assert_eq!(mmap.len(), 4);
}

#[test]
fn count_kmers_mmap_basic() {
    let path = fixture_path("simple.fa");
    let counts = count_kmers_mmap(&path, 4).expect("should count k-mers via mmap");

    assert!(!counts.is_empty(), "should find k-mers");
}

#[test]
fn count_kmers_mmap_matches_regular() {
    let path = fixture_path("simple.fa");

    let regular_counts = count_kmers(&path, 4).expect("regular count should work");
    let mmap_counts = count_kmers_mmap(&path, 4).expect("mmap count should work");

    assert_eq!(
        regular_counts.len(),
        mmap_counts.len(),
        "should have same number of k-mers"
    );

    for (kmer, count) in &regular_counts {
        assert_eq!(
            mmap_counts.get(kmer),
            Some(count),
            "k-mer '{kmer}' should have same count"
        );
    }
}

#[test]
fn count_kmers_mmap_handles_with_n() {
    let path = fixture_path("with_n.fa");
    let counts = count_kmers_mmap(&path, 4).expect("should handle sequences with N");

    // Should skip k-mers containing N but still find valid ones
    assert!(!counts.is_empty(), "should find valid k-mers");
}

#[test]
fn mmap_fasta_nonexistent_file_error() {
    let result = MmapFasta::open("nonexistent_file.fa");
    assert!(result.is_err(), "should error on nonexistent file");
}
