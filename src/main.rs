use std::process;

use colored::Colorize;
use krust::{cli, config::Config, run};

fn main() {
    let matches = cli::cli().get_matches();

    let k = matches.get_one::<String>("k").expect("required");
    let path = matches.get_one::<String>("path").expect("required");

    let config = Config::new(k, path).unwrap_or_else(|e| {
        println!();
        println!(
            "{}\n {}",
            "Problem parsing arguments:".blue().bold(),
            e.to_string().blue()
        );
        println!();
        println!(
            "{}\n {}\n  {}\n   {}",
            "Help menu:".blue().bold(),
            "$ cargo run -- --help".bold(),
            "or".underline(),
            "$ krust --help".bold()
        );
        println!();
        process::exit(1);
    });

    println!("{}: {}", "k-length".bold(), k.blue().bold());
    println!("{}: {}", "data".bold(), path.underline().bold().blue());
    println!(
        "{}: {}",
        "reader".bold(),
        match cfg!(feature = "needletail") {
            true => "needletail",
            _ => "rust-bio",
        }
        .blue()
        .bold()
    );
    println!();

    if let Err(e) = run::run(config.path, config.k) {
        eprintln!(
            "{}\n {}",
            "Application error:".blue().bold(),
            e.to_string().blue()
        );
        drop(e);
        process::exit(1);
    }
}
