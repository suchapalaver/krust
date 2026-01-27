#![allow(
    clippy::too_many_lines,
    clippy::needless_pass_by_value,
    clippy::expect_used,
    clippy::unwrap_used,
    clippy::redundant_clone
)]

use std::process;

use clap::Parser;
use colored::Colorize;
use kmerust::{
    cli::{Args, Cli, Command, QueryArgs},
    format::SequenceFormat,
    index::{load_index, save_index, KmerIndex},
    input::Input,
    kmer::{Kmer, KmerLength},
    run::{self, output_counts},
};

/// Initialize the tracing subscriber with environment filter.
///
/// Set `RUST_LOG=kmerust=debug` to see debug output.
#[cfg(feature = "tracing")]
fn init_tracing() {
    use tracing_subscriber::EnvFilter;

    tracing_subscriber::fmt()
        .with_env_filter(EnvFilter::from_default_env())
        .init();
}

fn main() {
    #[cfg(feature = "tracing")]
    init_tracing();

    // Check if first argument is a subcommand
    let args: Vec<String> = std::env::args().collect();
    if args.get(1).is_some_and(|arg| arg == "query") {
        // Parse as subcommand
        let cli = Cli::parse();
        if let Some(Command::Query(query_args)) = cli.command {
            run_query(query_args);
        }
        return;
    }

    // Otherwise, parse as count command (backward compatible)
    let args = Args::parse();
    run_count(args);
}

fn run_count(args: Args) {
    let input = args.input();

    // Validate that the file exists (only for file inputs)
    if let Input::File(ref path) = input {
        if !path.exists() {
            eprintln!(
                "{}\n {}",
                "Problem with arguments:".blue().bold(),
                format!("File not found: {}", path.display()).blue().bold()
            );
            process::exit(1);
        }
    }

    // Resolve input format (auto-detect from extension or use explicit setting)
    let input_format = match &input {
        Input::File(path) => args.input_format.resolve(Some(path)),
        Input::Stdin => args.input_format.resolve(None),
    };

    if !args.quiet {
        eprintln!(
            "{}: {}",
            "k-length".bold(),
            args.k.to_string().blue().bold()
        );
        eprintln!(
            "{}: {}",
            "data".bold(),
            input.to_string().underline().bold().blue()
        );
        // Show input format (detected or explicit)
        let format_display = match (args.input_format, input_format) {
            (SequenceFormat::Auto, resolved) => format!("{resolved} (auto-detected)"),
            (explicit, _) => explicit.to_string(),
        };
        eprintln!(
            "{}: {}",
            "input-format".bold(),
            format_display.blue().bold()
        );
        eprintln!(
            "{}: {}",
            "reader".bold(),
            if cfg!(feature = "needletail") {
                "needletail"
            } else {
                "rust-bio"
            }
            .blue()
            .bold()
        );
        eprintln!(
            "{}: {}",
            "output-format".bold(),
            format!("{:?}", args.format).to_lowercase().blue().bold()
        );
        if args.min_count > 1 {
            eprintln!(
                "{}: {}",
                "min-count".bold(),
                args.min_count.to_string().blue().bold()
            );
        }
        if let Some(min_q) = args.min_quality {
            eprintln!(
                "{}: {}",
                "min-quality".bold(),
                min_q.to_string().blue().bold()
            );
        }
        if let Some(ref save_path) = args.save {
            eprintln!(
                "{}: {}",
                "save-index".bold(),
                save_path.display().to_string().blue().bold()
            );
        }
        eprintln!();
    }

    // Warn if min-quality is set but input is FASTA (quality scores not available)
    if args.min_quality.is_some() && input_format.is_fasta() {
        eprintln!(
            "{}: {}",
            "warning".yellow().bold(),
            "--min-quality is ignored for FASTA input".yellow()
        );
    }

    // Warn if min-quality is set with stdin (not yet supported for stdin)
    if args.min_quality.is_some() && matches!(input, Input::Stdin) {
        eprintln!(
            "{}: {}",
            "warning".yellow().bold(),
            "--min-quality is not yet supported for stdin input".yellow()
        );
    }

    // If saving to index, we need to capture the packed counts
    if let Some(ref save_path) = args.save {
        // Count k-mers and get packed representation for the index
        let counts = match &input {
            Input::File(path) => {
                run::count_kmers_with_quality(path, args.k, input_format, args.min_quality)
                    .unwrap_or_else(|e| {
                        eprintln!(
                            "{}\n {}",
                            "Application error:".blue().bold(),
                            e.to_string().blue()
                        );
                        process::exit(1);
                    })
            }
            Input::Stdin => kmerust::streaming::count_kmers_stdin_with_format(args.k, input_format)
                .unwrap_or_else(|e| {
                    eprintln!(
                        "{}\n {}",
                        "Application error:".blue().bold(),
                        e.to_string().blue()
                    );
                    process::exit(1);
                }),
        };

        // Convert to packed counts for the index
        let k_len = KmerLength::new(args.k).unwrap();
        let packed_counts = counts_to_packed(&counts);

        // Save the index
        let index = KmerIndex::new(k_len, packed_counts);
        if let Err(e) = save_index(&index, save_path) {
            eprintln!(
                "{}\n {}",
                "Failed to save index:".blue().bold(),
                e.to_string().blue()
            );
            process::exit(1);
        }

        if !args.quiet {
            eprintln!(
                "{}: {} ({} k-mers)",
                "saved".bold(),
                save_path.display().to_string().green().bold(),
                index.len()
            );
        }

        // Also output to stdout as usual
        if let Err(e) = output_counts(counts, args.format, args.min_count) {
            eprintln!(
                "{}\n {}",
                "Application error:".blue().bold(),
                e.to_string().blue()
            );
            process::exit(1);
        }
    } else {
        // Normal operation: count and output
        if let Err(e) = run::run_with_quality(
            &input,
            args.k,
            args.format,
            args.min_count,
            input_format,
            args.min_quality,
        ) {
            eprintln!(
                "{}\n {}",
                "Application error:".blue().bold(),
                e.to_string().blue()
            );
            process::exit(1);
        }
    }
}

fn run_query(args: QueryArgs) {
    // Load the index
    let index = match load_index(&args.index) {
        Ok(idx) => idx,
        Err(e) => {
            eprintln!(
                "{}\n {}",
                "Failed to load index:".blue().bold(),
                e.to_string().blue()
            );
            process::exit(1);
        }
    };

    // Validate and pack the query k-mer
    let kmer_bytes = args.kmer.to_uppercase();
    if kmer_bytes.len() != index.k().get() {
        eprintln!(
            "{}\n {}",
            "Query error:".blue().bold(),
            format!(
                "k-mer length mismatch: query has {} bases, index has k={}",
                kmer_bytes.len(),
                index.k().get()
            )
            .blue()
        );
        process::exit(1);
    }

    // Pack the k-mer to its canonical form
    let packed = match Kmer::from_sub(bytes::Bytes::from(kmer_bytes.clone())) {
        Ok(kmer) => kmer.pack().canonical().packed_bits(),
        Err(e) => {
            eprintln!(
                "{}\n {}",
                "Invalid k-mer:".blue().bold(),
                e.to_string().blue()
            );
            process::exit(1);
        }
    };

    // Query and output
    match index.get(packed) {
        Some(count) => println!("{count}"),
        None => println!("0"),
    }
}

/// Convert string-keyed counts to packed representation.
fn counts_to_packed(
    counts: &std::collections::HashMap<String, u64>,
) -> std::collections::HashMap<u64, u64> {
    use bytes::Bytes;

    counts
        .iter()
        .map(|(kmer, &count)| {
            let packed = Kmer::from_sub(Bytes::from(kmer.clone()))
                .expect("k-mer from count should be valid")
                .pack()
                .packed_bits();
            (packed, count)
        })
        .collect()
}
