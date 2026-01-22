//! K-mer counting with progress reporting.
//!
//! This example demonstrates how to use the progress callback API to display
//! a progress indicator during long-running k-mer counting operations.
//!
//! # Usage
//!
//! ```bash
//! cargo run --example progress_bar -- large_genome.fa
//! ```

use std::env;
use std::io::{self, Write};
use std::process;
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::Arc;
use std::time::Instant;

use kmerust::builder::KmerCounter;

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() < 2 {
        eprintln!("Usage: {} <fasta_file> [k]", args[0]);
        eprintln!();
        eprintln!("Demonstrates progress reporting during k-mer counting.");
        process::exit(1);
    }

    let path = &args[1];
    let k: usize = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(21);

    eprintln!("Counting {k}-mers in {path}...\n");

    // Track the last reported values to avoid redundant updates
    let last_seqs = Arc::new(AtomicU64::new(0));
    let last_seqs_clone = Arc::clone(&last_seqs);

    let start = Instant::now();

    // Build counter and count with progress
    let counter = match KmerCounter::new().k(k) {
        Ok(c) => c,
        Err(e) => {
            eprintln!("Error: {e}");
            process::exit(1);
        }
    };

    let counts = match counter.count_with_progress(path, move |progress| {
        // Only update every 100 sequences to reduce overhead
        let prev = last_seqs_clone.load(Ordering::Relaxed);
        if progress.sequences_processed >= prev + 100 || progress.sequences_processed < prev {
            last_seqs_clone.store(progress.sequences_processed, Ordering::Relaxed);

            // Format bases with SI prefix
            let bases_str = format_bases(progress.bases_processed);

            // Print progress (overwrite same line)
            eprint!(
                "\r  Sequences: {:>8}  |  Bases: {:>10}",
                progress.sequences_processed, bases_str
            );
            let _ = io::stderr().flush();
        }
    }) {
        Ok(counts) => counts,
        Err(e) => {
            eprintln!("\nError: {e}");
            process::exit(1);
        }
    };

    let elapsed = start.elapsed();

    // Clear the progress line and print summary
    eprintln!("\r{:60}", ""); // Clear line
    eprintln!("\n=== Results ===");
    eprintln!("Unique k-mers:    {}", counts.len());
    eprintln!("Processing time:  {:.2?}", elapsed);

    // Calculate throughput
    let total_bases: u64 = last_seqs.load(Ordering::Relaxed);
    if elapsed.as_secs_f64() > 0.0 && total_bases > 0 {
        let bases_per_sec = total_bases as f64 / elapsed.as_secs_f64();
        eprintln!(
            "Throughput:       {} bases/sec",
            format_bases(bases_per_sec as u64)
        );
    }

    // Show top k-mers
    let mut sorted: Vec<_> = counts.into_iter().collect();
    sorted.sort_by(|a, b| b.1.cmp(&a.1));

    eprintln!("\nTop 10 k-mers:");
    for (kmer, count) in sorted.into_iter().take(10) {
        eprintln!("  {kmer}: {count}");
    }
}

/// Format a base count with SI prefix (K, M, G).
fn format_bases(bases: u64) -> String {
    if bases >= 1_000_000_000 {
        format!("{:.2}G", bases as f64 / 1_000_000_000.0)
    } else if bases >= 1_000_000 {
        format!("{:.2}M", bases as f64 / 1_000_000.0)
    } else if bases >= 1_000 {
        format!("{:.2}K", bases as f64 / 1_000.0)
    } else {
        format!("{bases}")
    }
}
