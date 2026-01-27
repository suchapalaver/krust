//! Basic k-mer counting example.
//!
//! This example demonstrates the simplest way to count k-mers using the builder API.
//!
//! # Usage
//!
//! ```bash
//! cargo run --example basic_count -- sequences.fa
//! ```

#![allow(clippy::unwrap_used, clippy::expect_used)]

use std::env;
use std::process;

use kmerust::builder::KmerCounter;

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() < 2 {
        eprintln!("Usage: {} <fasta_file> [k]", args[0]);
        eprintln!();
        eprintln!("Arguments:");
        eprintln!("  fasta_file  Path to a FASTA file");
        eprintln!("  k           K-mer length (default: 21)");
        process::exit(1);
    }

    let path = &args[1];
    let k: usize = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(21);

    // Count k-mers using the builder API
    let counts = match KmerCounter::new().k(k) {
        Ok(counter) => match counter.count(path) {
            Ok(counts) => counts,
            Err(e) => {
                eprintln!("Error counting k-mers: {e}");
                process::exit(1);
            }
        },
        Err(e) => {
            eprintln!("Invalid k-mer length: {e}");
            process::exit(1);
        }
    };

    // Print summary statistics
    println!("K-mer counting complete!");
    println!("  K-mer length: {k}");
    println!("  Unique k-mers: {}", counts.len());

    // Print the top 10 most frequent k-mers
    let mut sorted: Vec<_> = counts.into_iter().collect();
    sorted.sort_by(|a, b| b.1.cmp(&a.1));

    println!("\nTop 10 most frequent k-mers:");
    for (kmer, count) in sorted.into_iter().take(10) {
        println!("  {kmer}: {count}");
    }
}
