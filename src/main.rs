use std::process;

use krust::{config::Config, startup};

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

    println!();
    let k = matches.get_one::<String>("k").expect("required");
    println!("counting {}-mers", k);
    let path = matches.get_one::<String>("path").expect("required");
    println!("in {}", path);
    let reader = matches.get_one::<String>("reader").unwrap();
    println!("using {}", reader);
    println!();

    let config = Config::new(k, path, reader).unwrap_or_else(|err| {
        eprintln!("Problem parsing arguments: {}", err);
        process::exit(1);
    });

    if let Err(e) = startup::run(config.path, config.k, config.reader) {
        eprintln!("Application error: {}", e);
        drop(e);
        process::exit(1);
    }
}
