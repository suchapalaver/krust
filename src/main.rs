use std::{env, process};

fn main() {
    let config = krust::configuration::Config::new(env::args()).unwrap_or_else(|err| {
        eprintln!("Problem parsing arguments: {}", err);
        process::exit(1);
    });

    if let Err(e) = krust::startup::run(config.path, config.k) {
        eprintln!("Application error: {}", e);
        drop(e);
        process::exit(1);
    }
}
