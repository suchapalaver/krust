//! Jellyfish compatibility tests.
//!
//! These tests verify that kmerust produces the same output as Jellyfish,
//! the reference implementation for k-mer counting.
//!
//! Tests are marked with `#[ignore]` by default since they require:
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::ignore_without_reason
)]
//!
//! 1. Jellyfish to be installed (`brew install jellyfish` or `apt install jellyfish`)
//! 2. Additional test fixtures
//!
//! Run with: `cargo test --test jellyfish_compat -- --ignored`

use kmerust::run::count_kmers;
use std::collections::HashMap;
use std::io::Write;
use std::process::Command;
use tempfile::NamedTempFile;

/// Check if jellyfish is available on the system.
fn jellyfish_available() -> bool {
    Command::new("jellyfish")
        .arg("--version")
        .output()
        .map(|o| o.status.success())
        .unwrap_or(false)
}

/// Run jellyfish count on a file and return the k-mer counts.
fn run_jellyfish(path: &str, k: usize) -> Result<HashMap<String, u64>, String> {
    // Create temp file for jellyfish database
    let db_file = NamedTempFile::new().map_err(|e| e.to_string())?;

    // Run jellyfish count with canonical mode (-C)
    let count_output = Command::new("jellyfish")
        .args([
            "count",
            "-m",
            &k.to_string(),
            "-s",
            "10M",
            "-C", // canonical mode
            "-o",
            db_file.path().to_str().unwrap(),
            path,
        ])
        .output()
        .map_err(|e| format!("Failed to run jellyfish count: {e}"))?;

    if !count_output.status.success() {
        return Err(format!(
            "jellyfish count failed: {}",
            String::from_utf8_lossy(&count_output.stderr)
        ));
    }

    // Run jellyfish dump to get the k-mers
    let dump_output = Command::new("jellyfish")
        .args(["dump", "-c", db_file.path().to_str().unwrap()])
        .output()
        .map_err(|e| format!("Failed to run jellyfish dump: {e}"))?;

    if !dump_output.status.success() {
        return Err(format!(
            "jellyfish dump failed: {}",
            String::from_utf8_lossy(&dump_output.stderr)
        ));
    }

    // Parse jellyfish output (format: "KMER COUNT\n")
    let stdout = String::from_utf8_lossy(&dump_output.stdout);
    let mut counts = HashMap::new();

    for line in stdout.lines() {
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() == 2 {
            let kmer = parts[0].to_string();
            let count: u64 = parts[1].parse().map_err(|e| format!("Parse error: {e}"))?;
            counts.insert(kmer, count);
        }
    }

    Ok(counts)
}

/// Creates a temporary FASTA file with the given content.
fn temp_fasta(content: &str) -> NamedTempFile {
    let mut file = NamedTempFile::new().expect("Failed to create temp file");
    file.write_all(content.as_bytes())
        .expect("Failed to write temp file");
    file.flush().expect("Failed to flush temp file");
    file
}

#[test]
#[ignore]
fn matches_jellyfish_simple() {
    if !jellyfish_available() {
        eprintln!("Skipping test: jellyfish not installed");
        return;
    }

    let fasta = temp_fasta(">seq\nACGTACGTACGT\n");
    let path = fasta.path().to_str().unwrap();

    let jf_counts = run_jellyfish(path, 5).expect("jellyfish failed");
    let our_counts = count_kmers(path, 5).expect("kmerust failed");

    // Compare counts
    assert_eq!(
        jf_counts.len(),
        our_counts.len(),
        "Different number of k-mers: jellyfish={}, kmerust={}",
        jf_counts.len(),
        our_counts.len()
    );

    for (kmer, jf_count) in &jf_counts {
        let our_count = our_counts.get(kmer);
        assert_eq!(
            Some(jf_count),
            our_count,
            "Mismatch for k-mer {kmer}: jellyfish={jf_count}, kmerust={our_count:?}"
        );
    }
}

#[test]
#[ignore]
fn matches_jellyfish_with_n_bases() {
    if !jellyfish_available() {
        eprintln!("Skipping test: jellyfish not installed");
        return;
    }

    let fasta = temp_fasta(">seq1\nACGTNACGT\n>seq2\nGATTACA\n");
    let path = fasta.path().to_str().unwrap();

    let jf_counts = run_jellyfish(path, 3).expect("jellyfish failed");
    let our_counts = count_kmers(path, 3).expect("kmerust failed");

    assert_eq!(
        jf_counts.len(),
        our_counts.len(),
        "Different number of k-mers"
    );

    for (kmer, jf_count) in &jf_counts {
        let our_count = our_counts.get(kmer);
        assert_eq!(
            Some(jf_count),
            our_count,
            "Mismatch for k-mer {kmer}: jellyfish={jf_count}, kmerust={our_count:?}"
        );
    }
}

#[test]
#[ignore]
fn matches_jellyfish_soft_masked() {
    if !jellyfish_available() {
        eprintln!("Skipping test: jellyfish not installed");
        return;
    }

    // Test with soft-masked (lowercase) bases
    let fasta = temp_fasta(">seq\nacgtACGTacgt\n");
    let path = fasta.path().to_str().unwrap();

    let jf_counts = run_jellyfish(path, 5).expect("jellyfish failed");
    let our_counts = count_kmers(path, 5).expect("kmerust failed");

    assert_eq!(
        jf_counts.len(),
        our_counts.len(),
        "Different number of k-mers"
    );

    for (kmer, jf_count) in &jf_counts {
        let our_count = our_counts.get(kmer);
        assert_eq!(
            Some(jf_count),
            our_count,
            "Mismatch for k-mer {kmer}: jellyfish={jf_count}, kmerust={our_count:?}"
        );
    }
}

#[test]
#[ignore]
fn matches_jellyfish_simple_fixture() {
    if !jellyfish_available() {
        eprintln!("Skipping test: jellyfish not installed");
        return;
    }

    let path = "tests/fixtures/simple.fa";

    for k in [3, 5, 7] {
        let jf_counts = run_jellyfish(path, k).expect("jellyfish failed");
        let our_counts = count_kmers(path, k).expect("kmerust failed");

        assert_eq!(
            jf_counts.len(),
            our_counts.len(),
            "k={k}: Different number of k-mers: jellyfish={}, kmerust={}",
            jf_counts.len(),
            our_counts.len()
        );

        for (kmer, jf_count) in &jf_counts {
            let our_count = our_counts.get(kmer);
            assert_eq!(
                Some(jf_count),
                our_count,
                "k={k}: Mismatch for k-mer {kmer}: jellyfish={jf_count}, kmerust={our_count:?}"
            );
        }
    }
}

#[test]
#[ignore]
fn matches_jellyfish_with_n_fixture() {
    if !jellyfish_available() {
        eprintln!("Skipping test: jellyfish not installed");
        return;
    }

    let path = "tests/fixtures/with_n.fa";

    for k in [3, 5] {
        let jf_counts = run_jellyfish(path, k).expect("jellyfish failed");
        let our_counts = count_kmers(path, k).expect("kmerust failed");

        assert_eq!(
            jf_counts.len(),
            our_counts.len(),
            "k={k}: Different number of k-mers: jellyfish={}, kmerust={}",
            jf_counts.len(),
            our_counts.len()
        );

        for (kmer, jf_count) in &jf_counts {
            let our_count = our_counts.get(kmer);
            assert_eq!(
                Some(jf_count),
                our_count,
                "k={k}: Mismatch for k-mer {kmer}: jellyfish={jf_count}, kmerust={our_count:?}"
            );
        }
    }
}

#[test]
#[ignore]
fn matches_jellyfish_k_boundaries() {
    if !jellyfish_available() {
        eprintln!("Skipping test: jellyfish not installed");
        return;
    }

    // Test boundary k values
    let fasta = temp_fasta(">seq\nACGTACGTACGTACGTACGTACGTACGTACGTACGT\n");
    let path = fasta.path().to_str().unwrap();

    for k in [1, 2, 31, 32] {
        let jf_counts = run_jellyfish(path, k).expect("jellyfish failed");
        let our_counts = count_kmers(path, k).expect("kmerust failed");

        assert_eq!(
            jf_counts.len(),
            our_counts.len(),
            "k={k}: Different number of k-mers: jellyfish={}, kmerust={}",
            jf_counts.len(),
            our_counts.len()
        );

        for (kmer, jf_count) in &jf_counts {
            let our_count = our_counts.get(kmer);
            assert_eq!(
                Some(jf_count),
                our_count,
                "k={k}: Mismatch for k-mer {kmer}: jellyfish={jf_count}, kmerust={our_count:?}"
            );
        }
    }
}

#[test]
#[ignore]
fn matches_jellyfish_homopolymer() {
    if !jellyfish_available() {
        eprintln!("Skipping test: jellyfish not installed");
        return;
    }

    // Homopolymer runs are a common edge case
    let fasta = temp_fasta(">seq\nAAAAAAAAAAAAAAAA\n");
    let path = fasta.path().to_str().unwrap();

    for k in [1, 3, 5, 7] {
        let jf_counts = run_jellyfish(path, k).expect("jellyfish failed");
        let our_counts = count_kmers(path, k).expect("kmerust failed");

        assert_eq!(
            jf_counts.len(),
            our_counts.len(),
            "k={k}: Different number of k-mers"
        );

        for (kmer, jf_count) in &jf_counts {
            let our_count = our_counts.get(kmer);
            assert_eq!(
                Some(jf_count),
                our_count,
                "k={k}: Mismatch for k-mer {kmer}: jellyfish={jf_count}, kmerust={our_count:?}"
            );
        }
    }
}
