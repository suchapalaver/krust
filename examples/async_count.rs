//! Async k-mer counting example.
//!
//! This example demonstrates how to use the async API for k-mer counting,
//! which is useful when integrating with async runtimes like Tokio.
//!
//! # Usage
//!
//! ```bash
//! cargo run --example async_count --features async -- sequences.fa
//! ```
//!
#![allow(clippy::unwrap_used, clippy::expect_used)]
//! # Feature Flag
//!
//! Requires the `async` feature to be enabled.

#[cfg(feature = "async")]
use kmerust::async_api::AsyncKmerCounter;

#[cfg(feature = "async")]
#[tokio::main]
async fn main() {
    use std::env;
    use std::process;

    let args: Vec<String> = env::args().collect();

    if args.len() < 2 {
        eprintln!("Usage: {} <fasta_file> [k]", args[0]);
        eprintln!();
        eprintln!("Arguments:");
        eprintln!("  fasta_file  Path to a FASTA file");
        eprintln!("  k           K-mer length (default: 21)");
        process::exit(1);
    }

    let path = args[1].clone();
    let k: usize = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(21);

    eprintln!("Counting {k}-mers in {path} (async mode)...");

    // Use the async builder API
    let counter = match AsyncKmerCounter::new().k(k) {
        Ok(c) => c,
        Err(e) => {
            eprintln!("Invalid k-mer length: {e}");
            process::exit(1);
        }
    };

    let counts = match counter.count(path).await {
        Ok(counts) => counts,
        Err(e) => {
            eprintln!("Error counting k-mers: {e}");
            process::exit(1);
        }
    };

    eprintln!("Found {} unique k-mers", counts.len());

    // Print top 10
    let mut sorted: Vec<_> = counts.into_iter().collect();
    sorted.sort_by(|a, b| b.1.cmp(&a.1));

    println!("\nTop 10 most frequent k-mers:");
    for (kmer, count) in sorted.into_iter().take(10) {
        println!("{kmer}\t{count}");
    }
}

#[cfg(not(feature = "async"))]
fn main() {
    eprintln!("This example requires the 'async' feature.");
    eprintln!("Run with: cargo run --example async_count --features async -- <fasta_file>");
    std::process::exit(1);
}
