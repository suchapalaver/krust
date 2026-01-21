use std::process;

use clap::Parser;
use colored::Colorize;
use krust::{cli::Args, run};

fn main() {
    let args = Args::parse();

    // Validate that the file exists
    if !args.path.exists() {
        eprintln!(
            "{}\n {}",
            "Problem with arguments:".blue().bold(),
            format!("File not found: {}", args.path.display())
                .blue()
                .bold()
        );
        process::exit(1);
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
            args.path.display().to_string().underline().bold().blue()
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

    if let Err(e) = run::run_with_options(&args.path, args.k, args.format, args.min_count) {
        eprintln!(
            "{}\n {}",
            "Application error:".blue().bold(),
            e.to_string().blue()
        );
        process::exit(1);
    }
}
