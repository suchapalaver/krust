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
