use std::process;

use colored::Colorize;
use krust::{config::Config, run};

use clap::{Arg, Command};

fn main() {
    let matches = Command::new("krust")
        .version("1.0")
        .author("Joseph L. <jlivesey@gmail.com>")
        .about("krust: counts k-mers, written in rust")
        .arg(
            Arg::new("k")
                .help("provides k length, e.g. 5")
                .required(true),
        )
        .arg(
            Arg::new("path")
                .help("path to a FASTA file, e.g. /home/lisa/bio/cerevisiae.pan.fa")
                .required(true),
        )
        .arg(
            Arg::new("reader")
                .help("select *rust-bio* or *needletail* as FASTA reader")
                .required(false)
                .default_value("rust-bio"),
        )
        .get_matches();

    let k = matches.get_one::<String>("k").expect("required");
    let path = matches.get_one::<String>("path").expect("required");
    let reader = matches.get_one::<String>("reader").unwrap();
    
    println!();

    let config = Config::new(k, path, reader).unwrap_or_else(|e| {
        eprintln!("{}\n {}", "Problem parsing arguments:".red().bold(), e.to_string().red());
        eprintln!();
        eprintln!("{}\n {}\n  {}\n   {}", "Help menu:".green().bold(), "$ cargo run -- --help".bold(), "or".underline(), "$ krust --help".bold());
        eprintln!();
        process::exit(1);
    });

    println!("{}: {}", "k-length".bold(), k.yellow().bold());
    println!("{}: {}", "data".bold(), path.underline().bold().yellow());
    println!("{}: {}", "reader".bold(), reader.yellow().bold());
    println!();

    if let Err(e) = run::run(config.path, config.k, config.reader) {
        eprintln!("{}\n {}", "Application error:".red().bold(), e.to_string().red());
        drop(e);
        process::exit(1);
    }
}
