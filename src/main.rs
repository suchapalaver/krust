use std::process;

use clap::Parser;
use colored::Colorize;
use kmerust::{cli::Args, input::Input, run};

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
    let args = Args::parse();
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
            "format".bold(),
            format!("{:?}", args.format).to_lowercase().blue().bold()
        );
        if args.min_count > 1 {
            eprintln!(
                "{}: {}",
                "min-count".bold(),
                args.min_count.to_string().blue().bold()
            );
        }
        eprintln!();
    }

    if let Err(e) = run::run_with_input(&input, args.k, args.format, args.min_count) {
        eprintln!(
            "{}\n {}",
            "Application error:".blue().bold(),
            e.to_string().blue()
        );
        process::exit(1);
    }
}
