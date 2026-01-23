use std::process;

use clap::Parser;
use colored::Colorize;
use kmerust::{cli::Args, format::SequenceFormat, input::Input, run};

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
        eprintln!();
    }

    if let Err(e) =
        run::run_with_input_format(&input, args.k, args.format, args.min_count, input_format)
    {
        eprintln!(
            "{}\n {}",
            "Application error:".blue().bold(),
            e.to_string().blue()
        );
        process::exit(1);
    }
}
