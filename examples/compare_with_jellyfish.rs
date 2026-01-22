//! Compare kmerust output with Jellyfish.
//!
//! This example counts k-mers with kmerust and optionally compares with Jellyfish
//! output if Jellyfish is installed on the system.
//!
//! # Usage
//!
//! ```bash
//! cargo run --example compare_with_jellyfish -- sequences.fa
//! ```
//!
//! # Requirements
//!
//! For comparison, Jellyfish must be installed and available in PATH:
//! - macOS: `brew install jellyfish`
//! - Linux: `apt install jellyfish` or from source

use std::collections::HashMap;
use std::env;
use std::io::{BufRead, BufReader};
use std::process::{self, Command, Stdio};

use kmerust::builder::KmerCounter;

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() < 2 {
        eprintln!("Usage: {} <fasta_file> [k]", args[0]);
        eprintln!();
        eprintln!("Counts k-mers and compares with Jellyfish if available.");
        process::exit(1);
    }

    let path = &args[1];
    let k: usize = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(21);

    // Count with kmerust
    eprintln!("Counting {k}-mers with kmerust...");
    let kmerust_counts = match KmerCounter::new().k(k) {
        Ok(counter) => match counter.count(path) {
            Ok(counts) => counts,
            Err(e) => {
                eprintln!("Error: {e}");
                process::exit(1);
            }
        },
        Err(e) => {
            eprintln!("Error: {e}");
            process::exit(1);
        }
    };

    eprintln!(
        "kmerust found {} unique canonical k-mers",
        kmerust_counts.len()
    );

    // Try to run Jellyfish for comparison
    if let Some(jf_counts) = run_jellyfish(path, k) {
        compare_counts(&kmerust_counts, &jf_counts);
    } else {
        eprintln!("\nJellyfish not found. Install it to enable comparison:");
        eprintln!("  macOS:  brew install jellyfish");
        eprintln!("  Linux:  apt install jellyfish");
        eprintln!("\nkmerust results only:");
        print_summary(&kmerust_counts);
    }
}

fn run_jellyfish(path: &str, k: usize) -> Option<HashMap<String, u64>> {
    // Check if jellyfish is available
    let check = Command::new("which").arg("jellyfish").output().ok()?;

    if !check.status.success() {
        return None;
    }

    eprintln!("Counting {k}-mers with Jellyfish...");

    // Create a temporary file for the jellyfish database
    let db_path = format!("/tmp/jf_compare_{}.jf", std::process::id());

    // Run jellyfish count
    let count_status = Command::new("jellyfish")
        .args([
            "count",
            "-m",
            &k.to_string(),
            "-s",
            "100M",
            "-C", // canonical mode
            "-o",
            &db_path,
            path,
        ])
        .status()
        .ok()?;

    if !count_status.success() {
        eprintln!("Jellyfish count failed");
        return None;
    }

    // Run jellyfish dump to get the counts
    let dump_output = Command::new("jellyfish")
        .args(["dump", "-c", &db_path])
        .stdout(Stdio::piped())
        .spawn()
        .ok()?;

    let stdout = dump_output.stdout?;
    let reader = BufReader::new(stdout);

    let mut counts = HashMap::new();
    for line in reader.lines().map_while(Result::ok) {
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() == 2 {
            if let Ok(count) = parts[1].parse::<u64>() {
                counts.insert(parts[0].to_string(), count);
            }
        }
    }

    // Clean up
    let _ = std::fs::remove_file(&db_path);

    eprintln!("Jellyfish found {} unique canonical k-mers", counts.len());
    Some(counts)
}

fn compare_counts(kmerust: &HashMap<String, u64>, jellyfish: &HashMap<String, u64>) {
    let mut mismatches = 0;
    let mut kmerust_only = 0;
    let mut jellyfish_only = 0;

    // Check kmerust counts against jellyfish
    for (kmer, &kr_count) in kmerust {
        match jellyfish.get(kmer) {
            Some(&jf_count) if kr_count != jf_count => {
                if mismatches < 5 {
                    eprintln!("  Mismatch: {kmer} kmerust={kr_count} jellyfish={jf_count}");
                }
                mismatches += 1;
            }
            None => {
                kmerust_only += 1;
            }
            _ => {}
        }
    }

    // Check for k-mers only in jellyfish
    for kmer in jellyfish.keys() {
        if !kmerust.contains_key(kmer) {
            jellyfish_only += 1;
        }
    }

    println!("\n=== Comparison Results ===");
    println!("kmerust unique k-mers:   {}", kmerust.len());
    println!("Jellyfish unique k-mers: {}", jellyfish.len());
    println!();

    if mismatches == 0 && kmerust_only == 0 && jellyfish_only == 0 {
        println!("PERFECT MATCH! All k-mer counts are identical.");
    } else {
        println!("Differences:");
        println!("  Count mismatches:    {mismatches}");
        println!("  Only in kmerust:     {kmerust_only}");
        println!("  Only in Jellyfish:   {jellyfish_only}");
    }
}

fn print_summary(counts: &HashMap<String, u64>) {
    let total: u64 = counts.values().sum();
    let max_count = counts.values().max().copied().unwrap_or(0);

    println!("\n=== K-mer Summary ===");
    println!("Unique k-mers: {}", counts.len());
    println!("Total k-mers:  {total}");
    println!("Max count:     {max_count}");

    // Show top 5
    let mut sorted: Vec<_> = counts.iter().collect();
    sorted.sort_by(|a, b| b.1.cmp(a.1));

    println!("\nTop 5 k-mers:");
    for (kmer, count) in sorted.into_iter().take(5) {
        println!("  {kmer}: {count}");
    }
}
