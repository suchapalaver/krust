//! Streaming k-mer counting for large files.
//!
//! This example demonstrates memory-efficient k-mer counting for very large FASTA files
//! using the streaming API, which processes sequences one at a time.
//!
//! # Usage
//!
//! ```bash
//! cargo run --example streaming_large_file -- large_genome.fa
//! ```

#![allow(clippy::unwrap_used, clippy::expect_used)]

use std::env;
use std::process;

use kmerust::builder::KmerCounter;

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() < 2 {
        eprintln!("Usage: {} <fasta_file> [k] [min_count]", args[0]);
        eprintln!();
        eprintln!("Arguments:");
        eprintln!("  fasta_file  Path to a FASTA file");
        eprintln!("  k           K-mer length (default: 21)");
        eprintln!("  min_count   Minimum count threshold (default: 1)");
        eprintln!();
        eprintln!("This example uses streaming I/O for memory-efficient processing.");
        process::exit(1);
    }

    let path = &args[1];
    let k: usize = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(21);
    let min_count: u64 = args.get(3).and_then(|s| s.parse().ok()).unwrap_or(1);

    eprintln!("Counting {k}-mers in {path} (streaming mode)...");
    eprintln!("Minimum count threshold: {min_count}");

    // Build the counter with streaming mode
    let counter = match KmerCounter::new().k(k) {
        Ok(c) => c.min_count(min_count),
        Err(e) => {
            eprintln!("Invalid k-mer length: {e}");
            process::exit(1);
        }
    };

    // Use streaming API for memory efficiency
    let counts = match counter.count_streaming(path) {
        Ok(counts) => counts,
        Err(e) => {
            eprintln!("Error counting k-mers: {e}");
            process::exit(1);
        }
    };

    eprintln!(
        "Found {} unique k-mers with count >= {min_count}",
        counts.len()
    );

    // Output in TSV format
    for (kmer, count) in &counts {
        println!("{kmer}\t{count}");
    }
}
