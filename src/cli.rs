use clap::{Arg, Command};

pub fn cli() -> Command {
    Command::new("krust")
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
}
