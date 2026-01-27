#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::similar_names,
    clippy::too_many_lines,
    clippy::stable_sort_primitive
)]

use std::process::Command;

fn kmerust_cmd() -> Command {
    Command::new(env!("CARGO_BIN_EXE_kmerust"))
}

#[test]
fn cli_help_flag() {
    let output = kmerust_cmd()
        .arg("--help")
        .output()
        .expect("Failed to execute");
    assert!(output.status.success());
    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(stdout.contains("kmerust"));
    assert!(stdout.contains("k-mers"));
}

#[test]
fn cli_version_flag() {
    let output = kmerust_cmd()
        .arg("--version")
        .output()
        .expect("Failed to execute");
    assert!(output.status.success());
    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(stdout.contains(env!("CARGO_PKG_VERSION")));
}

#[test]
fn cli_missing_args() {
    let output = kmerust_cmd().output().expect("Failed to execute");
    assert!(!output.status.success());
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("required") || stderr.contains("Usage"));
}

#[test]
fn cli_stdin_default_when_path_omitted() {
    // When path is omitted, stdin is used as input
    // Provide a valid FASTA header with empty sequence
    use std::io::Write;
    use std::process::Stdio;

    let mut child = kmerust_cmd()
        .arg("5")
        .arg("--quiet")
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .expect("Failed to spawn");

    // Write a minimal valid FASTA (header only, no sequence long enough for k=5)
    child
        .stdin
        .as_mut()
        .unwrap()
        .write_all(b">seq\nACGT\n")
        .expect("Failed to write to stdin");

    let output = child.wait_with_output().expect("Failed to wait");
    // Should succeed with valid input (even if no k-mers are long enough)
    assert!(output.status.success());
    // Output should be empty since sequence is shorter than k
    assert!(output.stdout.is_empty());
}

#[test]
fn cli_invalid_k() {
    let output = kmerust_cmd()
        .args(["abc", "tests/fixtures/simple.fa"])
        .output()
        .expect("Failed to execute");
    assert!(!output.status.success());
}

#[test]
fn cli_k_zero() {
    let output = kmerust_cmd()
        .args(["0", "tests/fixtures/simple.fa"])
        .output()
        .expect("Failed to execute");
    assert!(!output.status.success());
}

#[test]
fn cli_k_too_large() {
    let output = kmerust_cmd()
        .args(["33", "tests/fixtures/simple.fa"])
        .output()
        .expect("Failed to execute");
    assert!(!output.status.success());
}

#[test]
fn cli_invalid_file_path() {
    let output = kmerust_cmd()
        .args(["5", "/nonexistent/path/to/file.fa"])
        .output()
        .expect("Failed to execute");
    assert!(!output.status.success());
}

#[test]
fn cli_simple_kmer_counting() {
    let output = kmerust_cmd()
        .args(["3", "tests/fixtures/simple.fa"])
        .output()
        .expect("Failed to execute");
    assert!(output.status.success());
    let stdout = String::from_utf8_lossy(&output.stdout);
    // Output should contain k-mers in FASTA-like format (>count\nkmer)
    assert!(stdout.contains('>'));
}

#[test]
fn cli_handles_n_bases() {
    let output = kmerust_cmd()
        .args(["3", "tests/fixtures/with_n.fa"])
        .output()
        .expect("Failed to execute");
    assert!(output.status.success());
    let stdout = String::from_utf8_lossy(&output.stdout);
    // Should still produce output (N bases are skipped)
    assert!(stdout.contains('>'));
    // Output should not contain N bases
    assert!(!stdout.contains('N'));
}

#[test]
fn cli_kmer_length_1() {
    let output = kmerust_cmd()
        .args(["1", "tests/fixtures/simple.fa"])
        .output()
        .expect("Failed to execute");
    assert!(output.status.success());
    let stdout = String::from_utf8_lossy(&output.stdout);
    // With k=1, we should only see single nucleotide k-mers (A, C, G, T)
    // and their canonical forms (A/T -> A, C/G -> C)
    assert!(stdout.contains('>'));
}

#[test]
fn cli_kmer_length_32() {
    let output = kmerust_cmd()
        .args(["32", "tests/fixtures/simple.fa"])
        .output()
        .expect("Failed to execute");
    // Should succeed even if there are no k-mers of length 32 in the test file
    assert!(output.status.success());
}

#[test]
fn cli_format_tsv() {
    let output = kmerust_cmd()
        .args(["3", "tests/fixtures/simple.fa", "--format", "tsv"])
        .output()
        .expect("Failed to execute");
    assert!(output.status.success());
    let stdout = String::from_utf8_lossy(&output.stdout);
    // TSV format should have tab-separated kmer and count
    assert!(stdout.contains('\t'));
    // Should not have FASTA-like format
    assert!(!stdout.contains('>'));
}

#[test]
fn cli_format_json() {
    let output = kmerust_cmd()
        .args(["3", "tests/fixtures/simple.fa", "--format", "json"])
        .output()
        .expect("Failed to execute");
    assert!(output.status.success());
    let stdout = String::from_utf8_lossy(&output.stdout);
    // JSON format should have array brackets and kmer/count keys
    assert!(stdout.contains('['));
    assert!(stdout.contains("kmer"));
    assert!(stdout.contains("count"));
}

#[test]
fn cli_min_count_filter() {
    // First count without filter to get baseline
    let output_no_filter = kmerust_cmd()
        .args([
            "3",
            "tests/fixtures/simple.fa",
            "--format",
            "tsv",
            "--quiet",
        ])
        .output()
        .expect("Failed to execute");
    let lines_no_filter = String::from_utf8_lossy(&output_no_filter.stdout)
        .lines()
        .count();

    // Now with a high min-count that should filter out most k-mers
    let output_filtered = kmerust_cmd()
        .args([
            "3",
            "tests/fixtures/simple.fa",
            "--format",
            "tsv",
            "--quiet",
            "--min-count",
            "1000",
        ])
        .output()
        .expect("Failed to execute");
    let lines_filtered = String::from_utf8_lossy(&output_filtered.stdout)
        .lines()
        .count();

    // Filtered output should have fewer lines
    assert!(lines_no_filter > 0, "Test fixture should produce k-mers");
    assert!(
        lines_filtered < lines_no_filter,
        "High min-count should filter out k-mers"
    );
}

#[test]
fn cli_quiet_flag() {
    let output_normal = kmerust_cmd()
        .args(["3", "tests/fixtures/simple.fa"])
        .output()
        .expect("Failed to execute");

    let output_quiet = kmerust_cmd()
        .args(["3", "tests/fixtures/simple.fa", "--quiet"])
        .output()
        .expect("Failed to execute");

    // Both should succeed
    assert!(output_normal.status.success());
    assert!(output_quiet.status.success());

    // Quiet mode should have no stderr output
    let stderr_quiet = String::from_utf8_lossy(&output_quiet.stderr);
    assert!(
        stderr_quiet.is_empty(),
        "Quiet mode should not produce stderr"
    );

    // Normal mode should have stderr output (info messages)
    let stderr_normal = String::from_utf8_lossy(&output_normal.stderr);
    assert!(
        !stderr_normal.is_empty(),
        "Normal mode should produce info on stderr"
    );
}

#[test]
fn cli_handles_soft_masked_bases() {
    // Test that soft-masked (lowercase) bases are counted like uppercase
    // The fixture has "AAAa" which should produce 2 counts of "AAA"
    let output = kmerust_cmd()
        .args([
            "3",
            "tests/fixtures/soft_masked.fa",
            "--format",
            "tsv",
            "--quiet",
        ])
        .output()
        .expect("Failed to execute");
    assert!(output.status.success());
    let stdout = String::from_utf8_lossy(&output.stdout);
    // Should contain AAA k-mer with count of 2
    assert!(stdout.contains("AAA\t2"));
}

#[test]
fn cli_stdin_with_fasta_content() {
    use std::io::Write;
    use std::process::Stdio;

    let fasta_content = b">seq1\nACGTACGT\n";

    let mut child = kmerust_cmd()
        .args(["4", "-", "--format", "tsv", "--quiet"])
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .expect("Failed to spawn");

    // Write FASTA content to stdin
    child
        .stdin
        .as_mut()
        .unwrap()
        .write_all(fasta_content)
        .expect("Failed to write to stdin");

    let output = child.wait_with_output().expect("Failed to wait");
    assert!(output.status.success());

    let stdout = String::from_utf8_lossy(&output.stdout);
    // Should contain k-mers from the sequence
    assert!(stdout.contains("ACGT") || stdout.contains("TACG") || stdout.contains("CGTA"));
}

#[test]
fn cli_stdin_pipe_simulation() {
    use std::io::Write;
    use std::process::Stdio;

    // Simulate: echo ">s\nACGT" | krust 2
    let fasta_content = b">s\nACGT\n";

    let mut child = kmerust_cmd()
        .args(["2", "--format", "tsv", "--quiet"])
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .expect("Failed to spawn");

    child
        .stdin
        .as_mut()
        .unwrap()
        .write_all(fasta_content)
        .expect("Failed to write to stdin");

    let output = child.wait_with_output().expect("Failed to wait");
    assert!(output.status.success());

    let stdout = String::from_utf8_lossy(&output.stdout);
    // With k=2 on "ACGT", we get: AC, CG, GT
    // Canonical forms: AC, CG, AC (GT's RC is AC)
    // So AC appears twice, CG once
    assert!(!stdout.is_empty());
}

#[test]
fn cli_stdin_json_output() {
    use std::io::Write;
    use std::process::Stdio;

    let fasta_content = b">seq\nAAAA\n";

    let mut child = kmerust_cmd()
        .args(["2", "-", "--format", "json", "--quiet"])
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .expect("Failed to spawn");

    child
        .stdin
        .as_mut()
        .unwrap()
        .write_all(fasta_content)
        .expect("Failed to write to stdin");

    let output = child.wait_with_output().expect("Failed to wait");
    assert!(output.status.success());

    let stdout = String::from_utf8_lossy(&output.stdout);
    // Should be valid JSON
    assert!(stdout.contains('['));
    assert!(stdout.contains("kmer"));
    assert!(stdout.contains("AA"));
}

#[test]
fn cli_stdin_multiple_sequences() {
    use std::io::Write;
    use std::process::Stdio;

    let fasta_content = b">seq1\nACGT\n>seq2\nTGCA\n";

    let mut child = kmerust_cmd()
        .args(["2", "-", "--format", "tsv", "--quiet"])
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .expect("Failed to spawn");

    child
        .stdin
        .as_mut()
        .unwrap()
        .write_all(fasta_content)
        .expect("Failed to write to stdin");

    let output = child.wait_with_output().expect("Failed to wait");
    assert!(output.status.success());

    let stdout = String::from_utf8_lossy(&output.stdout);
    // Should have processed both sequences
    assert!(!stdout.is_empty());
}

// FASTQ Support Tests

#[test]
fn cli_simple_fastq_counting() {
    let output = kmerust_cmd()
        .args(["3", "tests/fixtures/simple.fq"])
        .output()
        .expect("Failed to execute");
    assert!(output.status.success());
    let stdout = String::from_utf8_lossy(&output.stdout);
    // Output should contain k-mers in FASTA-like format (>count\nkmer)
    assert!(stdout.contains('>'));
}

#[test]
fn cli_fastq_handles_n_bases() {
    let output = kmerust_cmd()
        .args(["3", "tests/fixtures/with_n.fq"])
        .output()
        .expect("Failed to execute");
    assert!(output.status.success());
    let stdout = String::from_utf8_lossy(&output.stdout);
    // Should still produce output (N bases are skipped)
    assert!(stdout.contains('>'));
    // Output should not contain N bases
    assert!(!stdout.contains('N'));
}

#[test]
fn cli_fastq_explicit_format() {
    let output = kmerust_cmd()
        .args(["3", "tests/fixtures/simple.fq", "--input-format", "fastq"])
        .output()
        .expect("Failed to execute");
    assert!(output.status.success());
    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(stdout.contains('>'));
}

#[test]
fn cli_stdin_fastq_explicit_format() {
    use std::io::Write;
    use std::process::Stdio;

    let fastq_content = b"@seq1\nACGTACGT\n+\nIIIIIIII\n";

    let mut child = kmerust_cmd()
        .args([
            "4",
            "-",
            "--input-format",
            "fastq",
            "--format",
            "tsv",
            "--quiet",
        ])
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .expect("Failed to spawn");

    child
        .stdin
        .as_mut()
        .unwrap()
        .write_all(fastq_content)
        .expect("Failed to write to stdin");

    let output = child.wait_with_output().expect("Failed to wait");
    assert!(output.status.success());

    let stdout = String::from_utf8_lossy(&output.stdout);
    // Should contain k-mers from the sequence
    assert!(stdout.contains("ACGT") || stdout.contains("TACG") || stdout.contains("CGTA"));
}

#[test]
fn cli_fastq_matches_fasta_counts() {
    // Count k-mers from FASTA and FASTQ files with same sequences
    let fasta_output = kmerust_cmd()
        .args([
            "3",
            "tests/fixtures/simple.fa",
            "--format",
            "tsv",
            "--quiet",
        ])
        .output()
        .expect("Failed to execute FASTA");

    let fastq_output = kmerust_cmd()
        .args([
            "3",
            "tests/fixtures/simple.fq",
            "--format",
            "tsv",
            "--quiet",
        ])
        .output()
        .expect("Failed to execute FASTQ");

    assert!(fasta_output.status.success());
    assert!(fastq_output.status.success());

    // Parse outputs and compare counts
    let fasta_counts = parse_tsv_counts(&fasta_output.stdout);
    let fastq_counts = parse_tsv_counts(&fastq_output.stdout);

    // Same sequences should produce same k-mer counts
    assert_eq!(
        fasta_counts, fastq_counts,
        "FASTA and FASTQ should produce identical k-mer counts"
    );
}

#[test]
fn cli_fastq_format_auto_detection() {
    // Test that format is auto-detected from .fq extension
    let output = kmerust_cmd()
        .args([
            "3",
            "tests/fixtures/simple.fq",
            "--format",
            "tsv",
            "--quiet",
        ])
        .output()
        .expect("Failed to execute");
    assert!(output.status.success());

    let stdout = String::from_utf8_lossy(&output.stdout);
    // Should have k-mer output (format was correctly detected as FASTQ)
    assert!(!stdout.is_empty());
}

#[test]
fn cli_input_format_flag_short() {
    // Test short flag -i for input format
    let output = kmerust_cmd()
        .args(["3", "tests/fixtures/simple.fq", "-i", "fastq", "--quiet"])
        .output()
        .expect("Failed to execute");
    assert!(output.status.success());
}

#[test]
#[cfg_attr(not(feature = "gzip"), ignore)]
fn cli_gzip_fastq_counting() {
    // Test gzip-compressed FASTQ file
    let output = kmerust_cmd()
        .args([
            "3",
            "tests/fixtures/simple.fq.gz",
            "--format",
            "tsv",
            "--quiet",
        ])
        .output()
        .expect("Failed to execute");
    assert!(output.status.success());

    let stdout = String::from_utf8_lossy(&output.stdout);
    // Should have k-mer output
    assert!(!stdout.is_empty());

    // Verify it produces same counts as uncompressed FASTQ
    let uncompressed_output = kmerust_cmd()
        .args([
            "3",
            "tests/fixtures/simple.fq",
            "--format",
            "tsv",
            "--quiet",
        ])
        .output()
        .expect("Failed to execute");

    let gzip_counts = parse_tsv_counts(&output.stdout);
    let plain_counts = parse_tsv_counts(&uncompressed_output.stdout);

    assert_eq!(
        gzip_counts, plain_counts,
        "Gzipped FASTQ should produce identical counts to uncompressed FASTQ"
    );
}

fn parse_tsv_counts(output: &[u8]) -> std::collections::HashMap<String, u64> {
    String::from_utf8_lossy(output)
        .lines()
        .filter_map(|line| {
            let mut parts = line.split('\t');
            let kmer = parts.next()?;
            let count: u64 = parts.next()?.parse().ok()?;
            Some((kmer.to_string(), count))
        })
        .collect()
}

// Histogram Output Tests

#[test]
fn cli_format_histogram() {
    let output = kmerust_cmd()
        .args([
            "3",
            "tests/fixtures/simple.fa",
            "--format",
            "histogram",
            "--quiet",
        ])
        .output()
        .expect("Failed to execute");
    assert!(output.status.success());

    let stdout = String::from_utf8_lossy(&output.stdout);
    // Histogram format should have tab-separated count and frequency
    assert!(stdout.contains('\t'));
    // Should not have FASTA-like format
    assert!(!stdout.contains('>'));
    // Should not have k-mer strings
    assert!(!stdout.contains("ACG") && !stdout.contains("CGT"));
}

#[test]
fn cli_histogram_tsv_format() {
    let output = kmerust_cmd()
        .args([
            "3",
            "tests/fixtures/simple.fa",
            "--format",
            "histogram",
            "--quiet",
        ])
        .output()
        .expect("Failed to execute");
    assert!(output.status.success());

    let stdout = String::from_utf8_lossy(&output.stdout);
    // Each line should be count\tfrequency (two numbers separated by tab)
    for line in stdout.lines() {
        let parts: Vec<&str> = line.split('\t').collect();
        assert_eq!(parts.len(), 2, "Histogram line should have two columns");
        assert!(
            parts[0].parse::<u64>().is_ok(),
            "First column should be count"
        );
        assert!(
            parts[1].parse::<u64>().is_ok(),
            "Second column should be frequency"
        );
    }
}

#[test]
fn cli_histogram_sums_correctly() {
    // Count k-mers normally
    let tsv_output = kmerust_cmd()
        .args([
            "3",
            "tests/fixtures/simple.fa",
            "--format",
            "tsv",
            "--quiet",
        ])
        .output()
        .expect("Failed to execute");
    assert!(tsv_output.status.success());

    // Count histogram
    let hist_output = kmerust_cmd()
        .args([
            "3",
            "tests/fixtures/simple.fa",
            "--format",
            "histogram",
            "--quiet",
        ])
        .output()
        .expect("Failed to execute");
    assert!(hist_output.status.success());

    // Number of unique k-mers should equal sum of histogram frequencies
    let tsv_lines = String::from_utf8_lossy(&tsv_output.stdout).lines().count();

    let hist_freq_sum: u64 = String::from_utf8_lossy(&hist_output.stdout)
        .lines()
        .filter_map(|line| {
            let parts: Vec<&str> = line.split('\t').collect();
            parts.get(1)?.parse::<u64>().ok()
        })
        .sum();

    assert_eq!(
        tsv_lines as u64, hist_freq_sum,
        "Sum of histogram frequencies should equal number of unique k-mers"
    );
}

#[test]
fn cli_histogram_with_min_count() {
    // Histogram with min-count filter
    let output = kmerust_cmd()
        .args([
            "3",
            "tests/fixtures/simple.fa",
            "--format",
            "histogram",
            "--min-count",
            "2",
            "--quiet",
        ])
        .output()
        .expect("Failed to execute");
    assert!(output.status.success());

    let stdout = String::from_utf8_lossy(&output.stdout);
    // All counts in histogram should be >= 2
    for line in stdout.lines() {
        let parts: Vec<&str> = line.split('\t').collect();
        if let Ok(count) = parts[0].parse::<u64>() {
            assert!(count >= 2, "Histogram counts should be >= min_count");
        }
    }
}

#[test]
fn cli_histogram_sorted_by_count() {
    let output = kmerust_cmd()
        .args([
            "3",
            "tests/fixtures/simple.fa",
            "--format",
            "histogram",
            "--quiet",
        ])
        .output()
        .expect("Failed to execute");
    assert!(output.status.success());

    let stdout = String::from_utf8_lossy(&output.stdout);
    let counts: Vec<u64> = stdout
        .lines()
        .filter_map(|line| {
            let parts: Vec<&str> = line.split('\t').collect();
            parts.first()?.parse::<u64>().ok()
        })
        .collect();

    // BTreeMap should produce sorted output
    let mut sorted_counts = counts.clone();
    sorted_counts.sort();
    assert_eq!(
        counts, sorted_counts,
        "Histogram output should be sorted by count"
    );
}

#[test]
fn cli_histogram_stdin() {
    use std::io::Write;
    use std::process::Stdio;

    let fasta_content = b">seq1\nAAAAAAAA\n";

    let mut child = kmerust_cmd()
        .args(["3", "-", "--format", "histogram", "--quiet"])
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .expect("Failed to spawn");

    child
        .stdin
        .as_mut()
        .unwrap()
        .write_all(fasta_content)
        .expect("Failed to write to stdin");

    let output = child.wait_with_output().expect("Failed to wait");
    assert!(output.status.success());

    let stdout = String::from_utf8_lossy(&output.stdout);
    // "AAAAAAAA" has 6 "AAA" k-mers, all the same (canonical)
    // So histogram should show count 6 with frequency 1
    assert!(
        stdout.contains("6\t1"),
        "Should have one k-mer with count 6"
    );
}

// ============================================================================
// Index Serialization Tests (--save flag and query command)
// ============================================================================

#[test]
fn cli_save_flag_creates_index() {
    let tmp = tempfile::NamedTempFile::with_suffix(".kmix").unwrap();

    let output = kmerust_cmd()
        .args([
            "4",
            "tests/fixtures/simple.fa",
            "--save",
            tmp.path().to_str().unwrap(),
            "--quiet",
        ])
        .output()
        .expect("Failed to execute");

    assert!(output.status.success(), "Command should succeed");
    assert!(tmp.path().exists(), "Index file should be created");

    // Verify the file is a valid index by checking magic bytes
    let data = std::fs::read(tmp.path()).unwrap();
    assert!(data.len() >= 18, "Index file should be at least 18 bytes");
    assert_eq!(&data[0..4], b"KMIX", "Index should have KMIX magic bytes");
}

#[test]
fn cli_save_flag_outputs_counts_and_saves() {
    let tmp = tempfile::NamedTempFile::with_suffix(".kmix").unwrap();

    let output = kmerust_cmd()
        .args([
            "4",
            "tests/fixtures/simple.fa",
            "--save",
            tmp.path().to_str().unwrap(),
            "--quiet",
            "--format",
            "tsv",
        ])
        .output()
        .expect("Failed to execute");

    assert!(output.status.success());
    assert!(tmp.path().exists());

    // Verify stdout has k-mer counts
    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(stdout.contains('\t'), "TSV output should contain tabs");
    assert!(!stdout.is_empty(), "Should have output k-mers");
}

#[test]
fn cli_query_basic() {
    let tmp = tempfile::NamedTempFile::with_suffix(".kmix").unwrap();

    // First create an index from simple.fa
    let create = kmerust_cmd()
        .args([
            "4",
            "tests/fixtures/simple.fa",
            "--save",
            tmp.path().to_str().unwrap(),
            "--quiet",
        ])
        .output()
        .expect("Failed to create index");
    assert!(create.status.success());

    // Now query a k-mer that should exist
    // simple.fa contains "ACGT" which has canonical k-mers
    let query = kmerust_cmd()
        .args(["query", tmp.path().to_str().unwrap(), "ACGT"])
        .output()
        .expect("Failed to query");

    assert!(query.status.success(), "Query should succeed");
    let count: u64 = String::from_utf8_lossy(&query.stdout)
        .trim()
        .parse()
        .expect("Output should be a number");
    assert!(count > 0, "ACGT should have a count > 0");
}

#[test]
fn cli_query_nonexistent_kmer() {
    let tmp = tempfile::NamedTempFile::with_suffix(".kmix").unwrap();

    // Create an index from simple.fa (small file)
    let create = kmerust_cmd()
        .args([
            "4",
            "tests/fixtures/simple.fa",
            "--save",
            tmp.path().to_str().unwrap(),
            "--quiet",
        ])
        .output()
        .expect("Failed to create index");
    assert!(create.status.success());

    // Query a k-mer that likely doesn't exist
    let query = kmerust_cmd()
        .args(["query", tmp.path().to_str().unwrap(), "ZZZZ"])
        .output()
        .expect("Failed to query");

    // Should fail because Z is not a valid base
    assert!(
        !query.status.success(),
        "Query with invalid base should fail"
    );
}

#[test]
fn cli_query_wrong_length() {
    let tmp = tempfile::NamedTempFile::with_suffix(".kmix").unwrap();

    // Create a k=4 index
    let create = kmerust_cmd()
        .args([
            "4",
            "tests/fixtures/simple.fa",
            "--save",
            tmp.path().to_str().unwrap(),
            "--quiet",
        ])
        .output()
        .expect("Failed to create index");
    assert!(create.status.success());

    // Query with wrong k-mer length (5 instead of 4)
    let query = kmerust_cmd()
        .args(["query", tmp.path().to_str().unwrap(), "ACGTA"])
        .output()
        .expect("Failed to query");

    assert!(
        !query.status.success(),
        "Query with wrong k-mer length should fail"
    );
    let stderr = String::from_utf8_lossy(&query.stderr);
    assert!(
        stderr.contains("mismatch") || stderr.contains("length"),
        "Should report length mismatch"
    );
}

#[test]
fn cli_query_case_insensitive() {
    let tmp = tempfile::NamedTempFile::with_suffix(".kmix").unwrap();

    // Create an index
    let create = kmerust_cmd()
        .args([
            "4",
            "tests/fixtures/simple.fa",
            "--save",
            tmp.path().to_str().unwrap(),
            "--quiet",
        ])
        .output()
        .expect("Failed to create index");
    assert!(create.status.success());

    // Query with lowercase
    let query_lower = kmerust_cmd()
        .args(["query", tmp.path().to_str().unwrap(), "acgt"])
        .output()
        .expect("Failed to query");

    // Query with uppercase
    let query_upper = kmerust_cmd()
        .args(["query", tmp.path().to_str().unwrap(), "ACGT"])
        .output()
        .expect("Failed to query");

    assert!(query_lower.status.success());
    assert!(query_upper.status.success());

    // Both should return the same count
    let count_lower: u64 = String::from_utf8_lossy(&query_lower.stdout)
        .trim()
        .parse()
        .unwrap();
    let count_upper: u64 = String::from_utf8_lossy(&query_upper.stdout)
        .trim()
        .parse()
        .unwrap();
    assert_eq!(count_lower, count_upper, "Case should not matter");
}

#[test]
fn cli_query_canonical_equivalence() {
    let tmp = tempfile::NamedTempFile::with_suffix(".kmix").unwrap();

    // Create an index
    let create = kmerust_cmd()
        .args([
            "4",
            "tests/fixtures/simple.fa",
            "--save",
            tmp.path().to_str().unwrap(),
            "--quiet",
        ])
        .output()
        .expect("Failed to create index");
    assert!(create.status.success());

    // Query a k-mer and its reverse complement - should give same count
    // ACGT reverse complement is ACGT (palindrome)
    // Let's try AAAA vs TTTT
    let query_a = kmerust_cmd()
        .args(["query", tmp.path().to_str().unwrap(), "AAAA"])
        .output()
        .expect("Failed to query");

    let query_t = kmerust_cmd()
        .args(["query", tmp.path().to_str().unwrap(), "TTTT"])
        .output()
        .expect("Failed to query");

    assert!(query_a.status.success());
    assert!(query_t.status.success());

    // AAAA and TTTT are reverse complements, so should have same count
    let count_a: u64 = String::from_utf8_lossy(&query_a.stdout)
        .trim()
        .parse()
        .unwrap();
    let count_t: u64 = String::from_utf8_lossy(&query_t.stdout)
        .trim()
        .parse()
        .unwrap();
    assert_eq!(
        count_a, count_t,
        "Reverse complements should have same count"
    );
}

#[test]
fn cli_query_invalid_index() {
    let tmp = tempfile::NamedTempFile::with_suffix(".kmix").unwrap();

    // Write garbage to the file
    std::fs::write(tmp.path(), b"NOT_A_VALID_INDEX_FILE").unwrap();

    let query = kmerust_cmd()
        .args(["query", tmp.path().to_str().unwrap(), "ACGT"])
        .output()
        .expect("Failed to query");

    assert!(
        !query.status.success(),
        "Query on invalid index should fail"
    );
    let stderr = String::from_utf8_lossy(&query.stderr);
    assert!(
        stderr.contains("invalid") || stderr.contains("Failed"),
        "Should report invalid index"
    );
}

#[test]
fn cli_query_help() {
    let output = kmerust_cmd()
        .args(["query", "--help"])
        .output()
        .expect("Failed to execute");

    assert!(output.status.success());
    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(
        stdout.contains("index") || stdout.contains("INDEX"),
        "Help should mention index"
    );
    assert!(
        stdout.contains("kmer") || stdout.contains("KMER"),
        "Help should mention kmer"
    );
}
