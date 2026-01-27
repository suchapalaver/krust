//! Direct library API tests.
//!
//! These tests call the library functions directly without going through the CLI,
//! enabling more precise assertions about behavior and return values.

#![allow(clippy::unwrap_used, clippy::expect_used)]

use kmerust::run::count_kmers;
use kmerust::streaming::count_kmers_from_reader;
use std::io::{BufReader, Write};
use tempfile::NamedTempFile;

/// Creates a temporary FASTA file with the given content and returns its path.
fn temp_fasta(content: &str) -> NamedTempFile {
    let mut file = NamedTempFile::new().expect("Failed to create temp file");
    file.write_all(content.as_bytes())
        .expect("Failed to write temp file");
    file.flush().expect("Failed to flush temp file");
    file
}

#[test]
fn count_kmers_basic() {
    let fasta = temp_fasta(">seq\nACGT\n");
    let result = count_kmers(fasta.path(), 3).unwrap();

    // ACGT has 2 3-mers: ACG, CGT
    // ACG -> reverse complement is CGT, ACG < CGT, so canonical is ACG
    // CGT -> reverse complement is ACG, ACG < CGT, so canonical is ACG
    // Both map to ACG!
    assert_eq!(result.get("ACG"), Some(&2));
    assert_eq!(result.len(), 1);
}

#[test]
fn count_kmers_simple_fixture() {
    // simple.fa contains:
    // >seq1
    // ACGTACGT
    // >seq2
    // GATTACA
    let result = count_kmers("tests/fixtures/simple.fa", 3).unwrap();

    // Verify we get results
    assert!(!result.is_empty());

    // All k-mers should have positive counts
    for (kmer, count) in &result {
        assert!(*count > 0, "k-mer {kmer} has non-positive count {count}");
        assert_eq!(kmer.len(), 3, "k-mer {kmer} is not length 3");
    }
}

#[test]
fn count_kmers_returns_canonical_forms() {
    // Create a sequence that has a k-mer and its reverse complement
    // TTT -> AAA (reverse complement), AAA is canonical
    let fasta = temp_fasta(">seq\nTTT\n");
    let result = count_kmers(fasta.path(), 3).unwrap();

    // Should return AAA (canonical form), not TTT
    assert_eq!(result.get("AAA"), Some(&1));
    assert_eq!(result.get("TTT"), None);
}

#[test]
fn count_kmers_handles_n_bases() {
    // N bases should cause k-mer windows to be skipped
    let fasta = temp_fasta(">seq\nACGNACG\n");
    let result = count_kmers(fasta.path(), 3).unwrap();

    // ACG appears twice (before and after N)
    // The N causes a gap, so we get ACG, CGN (skip), GNA (skip), NAC (skip), ACG
    assert_eq!(result.get("ACG"), Some(&2));

    // No k-mers should contain N
    for kmer in result.keys() {
        assert!(!kmer.contains('N'), "k-mer {kmer} should not contain N");
    }
}

#[test]
fn count_kmers_soft_masked_bases() {
    // Lowercase bases should be treated same as uppercase
    let fasta = temp_fasta(">seq\nacgt\n");
    let result = count_kmers(fasta.path(), 3).unwrap();

    // Same as ACGT - should get ACG with count 2
    assert_eq!(result.get("ACG"), Some(&2));
}

#[test]
fn count_kmers_mixed_case() {
    // Mixed case should work
    let fasta = temp_fasta(">seq\nAcGt\n");
    let result = count_kmers(fasta.path(), 3).unwrap();

    assert_eq!(result.get("ACG"), Some(&2));
}

#[test]
fn count_kmers_empty_result_for_short_sequence() {
    // Sequence shorter than k should produce no k-mers
    let fasta = temp_fasta(">seq\nAC\n");
    let result = count_kmers(fasta.path(), 3).unwrap();

    assert!(result.is_empty());
}

#[test]
fn count_kmers_exact_length_sequence() {
    // Sequence exactly k length should produce one k-mer
    let fasta = temp_fasta(">seq\nACG\n");
    let result = count_kmers(fasta.path(), 3).unwrap();

    assert_eq!(result.len(), 1);
    assert_eq!(result.get("ACG"), Some(&1));
}

#[test]
fn count_kmers_multiple_sequences() {
    // Multiple sequences should be combined
    let fasta = temp_fasta(">seq1\nACG\n>seq2\nACG\n");
    let result = count_kmers(fasta.path(), 3).unwrap();

    assert_eq!(result.get("ACG"), Some(&2));
}

#[test]
fn count_kmers_k_equals_1() {
    // k=1 should work
    let fasta = temp_fasta(">seq\nACGT\n");
    let result = count_kmers(fasta.path(), 1).unwrap();

    // A/T are complements, C/G are complements
    // A and T both map to A (canonical)
    // C and G both map to C (canonical)
    assert_eq!(result.get("A"), Some(&2)); // A + T
    assert_eq!(result.get("C"), Some(&2)); // C + G
}

#[test]
fn count_kmers_k_equals_32() {
    // k=32 (max) should work
    let seq = "A".repeat(32);
    let fasta = temp_fasta(&format!(">seq\n{seq}\n"));
    let result = count_kmers(fasta.path(), 32).unwrap();

    // AAAA...A (32 As) has canonical form AAAA...A
    assert_eq!(result.len(), 1);
    assert!(result.contains_key(&seq));
}

#[test]
fn count_kmers_rejects_k_zero() {
    let fasta = temp_fasta(">seq\nACGT\n");
    let result = count_kmers(fasta.path(), 0);

    assert!(result.is_err());
}

#[test]
fn count_kmers_rejects_k_too_large() {
    let fasta = temp_fasta(">seq\nACGT\n");
    let result = count_kmers(fasta.path(), 33);

    assert!(result.is_err());
}

#[test]
fn count_kmers_nonexistent_file() {
    let result = count_kmers("/nonexistent/path/to/file.fa", 3);

    assert!(result.is_err());
}

// needletail returns an error for empty files, while rust-bio returns empty results
#[test]
#[cfg(not(feature = "needletail"))]
fn count_kmers_empty_file() {
    let fasta = temp_fasta("");
    let result = count_kmers(fasta.path(), 3).unwrap();

    assert!(result.is_empty());
}

// needletail returns an error for header-only files, while rust-bio returns empty results
#[test]
#[cfg(not(feature = "needletail"))]
fn count_kmers_header_only() {
    // File with only a header and no sequence
    let fasta = temp_fasta(">seq\n");
    let result = count_kmers(fasta.path(), 3).unwrap();

    assert!(result.is_empty());
}

#[test]
fn count_kmers_palindrome() {
    // ACGT is its own reverse complement
    let fasta = temp_fasta(">seq\nACGT\n");
    let result = count_kmers(fasta.path(), 4).unwrap();

    // ACGT is palindromic, so canonical form is ACGT
    assert_eq!(result.get("ACGT"), Some(&1));
}

#[test]
fn count_kmers_all_same_base() {
    // Homopolymer run
    let fasta = temp_fasta(">seq\nAAAAA\n");
    let result = count_kmers(fasta.path(), 3).unwrap();

    // AAA -> TTT reverse complement, AAA < TTT, so AAA is canonical
    // AAAAA has 3 3-mers: AAA, AAA, AAA
    assert_eq!(result.get("AAA"), Some(&3));
}

#[test]
fn count_kmers_complementary_sequences() {
    // Sequence and its reverse complement should give same canonical k-mers
    let fasta1 = temp_fasta(">seq\nGATTACA\n");
    let fasta2 = temp_fasta(">seq\nTGTAATC\n"); // reverse complement of GATTACA

    let result1 = count_kmers(fasta1.path(), 3).unwrap();
    let result2 = count_kmers(fasta2.path(), 3).unwrap();

    // Canonical k-mers should be the same
    assert_eq!(result1, result2);
}

#[test]
fn count_kmers_multiline_sequence() {
    // FASTA with sequence split across multiple lines
    let fasta = temp_fasta(">seq\nACG\nTAC\n");
    let result = count_kmers(fasta.path(), 3).unwrap();

    // Should treat as ACGTAC (6 bases, 4 3-mers)
    // ACG, CGT, GTA, TAC
    assert!(!result.is_empty());
}

#[test]
fn count_kmers_with_n_fixture() {
    // with_n.fa contains:
    // >seq1
    // ACGTNACGT
    // >seq2
    // NNNGATTACANNN
    let result = count_kmers("tests/fixtures/with_n.fa", 3).unwrap();

    // Should have results from valid portions
    assert!(!result.is_empty());

    // No k-mers should contain N
    for kmer in result.keys() {
        assert!(!kmer.contains('N'), "k-mer should not contain N");
    }
}

#[test]
fn count_kmers_soft_masked_fixture() {
    // soft_masked.fa contains:
    // >seq
    // AAAa
    let result = count_kmers("tests/fixtures/soft_masked.fa", 3).unwrap();

    // AAAa treated as AAAA, has 2 3-mers: AAA, AAA
    assert_eq!(result.get("AAA"), Some(&2));
}

// Tests for count_kmers_from_reader (stdin/reader support)

#[test]
fn count_kmers_from_reader_basic() {
    let fasta_data = b">seq\nACGT\n";
    let reader = BufReader::new(&fasta_data[..]);
    let result = count_kmers_from_reader(reader, 3).unwrap();

    // Same as file-based test
    assert_eq!(result.get("ACG"), Some(&2));
    assert_eq!(result.len(), 1);
}

#[test]
fn count_kmers_from_reader_empty_sequence() {
    // Test with a header but no sequence (or sequence too short for k)
    let fasta_data = b">seq\nAC\n";
    let reader = BufReader::new(&fasta_data[..]);
    let result = count_kmers_from_reader(reader, 3).unwrap();

    // Sequence "AC" is too short for k=3, so no k-mers
    assert!(result.is_empty());
}

#[test]
fn count_kmers_from_reader_multiple_sequences() {
    let fasta_data = b">seq1\nACGT\n>seq2\nTGCA\n";
    let reader = BufReader::new(&fasta_data[..]);
    let result = count_kmers_from_reader(reader, 2).unwrap();

    // Should process both sequences
    assert!(!result.is_empty());
}

#[test]
fn count_kmers_from_reader_handles_n_bases() {
    let fasta_data = b">seq\nACGNACG\n";
    let reader = BufReader::new(&fasta_data[..]);
    let result = count_kmers_from_reader(reader, 3).unwrap();

    // ACG appears twice (before and after N)
    assert_eq!(result.get("ACG"), Some(&2));

    // No k-mers should contain N
    for kmer in result.keys() {
        assert!(!kmer.contains('N'), "k-mer {kmer} should not contain N");
    }
}

#[test]
fn count_kmers_from_reader_soft_masked() {
    let fasta_data = b">seq\nacgt\n";
    let reader = BufReader::new(&fasta_data[..]);
    let result = count_kmers_from_reader(reader, 3).unwrap();

    assert_eq!(result.get("ACG"), Some(&2));
}

#[test]
fn count_kmers_from_reader_matches_file_based() {
    // Verify that reader-based counting gives same results as file-based
    let fasta_content = ">seq1\nACGTACGT\n>seq2\nGATTACA\n";

    // File-based
    let fasta_file = temp_fasta(fasta_content);
    let file_result = count_kmers(fasta_file.path(), 3).unwrap();

    // Reader-based
    let reader = BufReader::new(fasta_content.as_bytes());
    let reader_result = count_kmers_from_reader(reader, 3).unwrap();

    assert_eq!(file_result, reader_result);
}

#[test]
fn count_kmers_from_reader_rejects_invalid_k() {
    let fasta_data = b">seq\nACGT\n";

    let reader = BufReader::new(&fasta_data[..]);
    assert!(count_kmers_from_reader(reader, 0).is_err());

    let reader = BufReader::new(&fasta_data[..]);
    assert!(count_kmers_from_reader(reader, 33).is_err());
}
