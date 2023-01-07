use std::{error::Error, fmt::Debug, path::Path};

use bio::io::fasta::Reader;
use bytes::Bytes;
use needletail::parse_fastx_file;
use rayon::{prelude::IntoParallelIterator, vec::IntoIter};

pub(crate) trait SequenceReader {
    fn sequence_reader<P: AsRef<Path> + Debug>(path: P) -> Result<IntoIter<Bytes>, Box<dyn Error>>;
}

pub(crate) struct RustBio;

impl SequenceReader for RustBio {
    fn sequence_reader<P: AsRef<Path> + Debug>(path: P) -> Result<IntoIter<Bytes>, Box<dyn Error>> {
        Ok(Reader::from_file(path)?
            .records()
            .into_iter()
            .map(|read| read.expect("Error reading fasta record."))
            .map(|record| Bytes::copy_from_slice(record.seq()))
            .collect::<Vec<Bytes>>()
            .into_par_iter())
    }
}

pub(crate) struct Needletail;

impl SequenceReader for Needletail {
    fn sequence_reader<P: AsRef<Path> + Debug>(path: P) -> Result<IntoIter<Bytes>, Box<dyn Error>> {
        let mut reader = parse_fastx_file(path)?;
        let mut v = Vec::new();
        while let Some(record) = reader.next() {
            let record = record.expect("invalid record");
            let seq = record.seq();
            let seq = Bytes::copy_from_slice(&seq);
            v.push(seq);
        }
        Ok(v.into_par_iter())
    }
}
