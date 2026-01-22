use std::process::Command;

fn krust_cmd() -> Command {
    Command::new(env!("CARGO_BIN_EXE_krust"))
}

#[test]
fn cli_help_flag() {
    let output = krust_cmd()
        .arg("--help")
        .output()
        .expect("Failed to execute");
    assert!(output.status.success());
    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(stdout.contains("krust"));
    assert!(stdout.contains("k-mers"));
}

#[test]
fn cli_version_flag() {
    let output = krust_cmd()
        .arg("--version")
        .output()
        .expect("Failed to execute");
    assert!(output.status.success());
    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(stdout.contains(env!("CARGO_PKG_VERSION")));
}

#[test]
fn cli_missing_args() {
    let output = krust_cmd().output().expect("Failed to execute");
    assert!(!output.status.success());
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("required") || stderr.contains("Usage"));
}

#[test]
fn cli_missing_path_arg() {
    let output = krust_cmd().arg("5").output().expect("Failed to execute");
    assert!(!output.status.success());
}

#[test]
fn cli_invalid_k() {
    let output = krust_cmd()
        .args(["abc", "tests/fixtures/simple.fa"])
        .output()
        .expect("Failed to execute");
    assert!(!output.status.success());
}

#[test]
fn cli_k_zero() {
    let output = krust_cmd()
        .args(["0", "tests/fixtures/simple.fa"])
        .output()
        .expect("Failed to execute");
    assert!(!output.status.success());
}

#[test]
fn cli_k_too_large() {
    let output = krust_cmd()
        .args(["33", "tests/fixtures/simple.fa"])
        .output()
        .expect("Failed to execute");
    assert!(!output.status.success());
}

#[test]
fn cli_invalid_file_path() {
    let output = krust_cmd()
        .args(["5", "/nonexistent/path/to/file.fa"])
        .output()
        .expect("Failed to execute");
    assert!(!output.status.success());
}

#[test]
fn cli_simple_kmer_counting() {
    let output = krust_cmd()
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
    let output = krust_cmd()
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
    let output = krust_cmd()
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
    let output = krust_cmd()
        .args(["32", "tests/fixtures/simple.fa"])
        .output()
        .expect("Failed to execute");
    // Should succeed even if there are no k-mers of length 32 in the test file
    assert!(output.status.success());
}

#[test]
fn cli_format_tsv() {
    let output = krust_cmd()
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
    let output = krust_cmd()
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
    let output_no_filter = krust_cmd()
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
    let output_filtered = krust_cmd()
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
    let output_normal = krust_cmd()
        .args(["3", "tests/fixtures/simple.fa"])
        .output()
        .expect("Failed to execute");

    let output_quiet = krust_cmd()
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
