use std::{env, process};

use krust::{config::Config, startup};

fn main() {
    let config = Config::new(env::args()).unwrap_or_else(|err| {
        eprintln!("Problem parsing arguments: {}", err);
        process::exit(1);
    });

    if let Err(e) = startup::run(config.path, config.k) {
        eprintln!("Application error: {}", e);
        drop(e);
        process::exit(1);
    }
}
