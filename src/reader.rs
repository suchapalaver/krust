use std::{error::Error, fmt::Debug, path::Path, vec::IntoIter};

use bio::io::fasta::Reader;
use bytes::Bytes;
use rayon::prelude::{ParallelBridge, ParallelIterator};

pub(crate) trait SequenceReader {
    fn sequence_reader<P: AsRef<Path> + Debug>(
        path: P,
    ) -> Result<IntoIter<Bytes>, Box<dyn Error>>;
}

pub(crate) struct RustBio;

impl SequenceReader for RustBio {
    fn sequence_reader<P: AsRef<Path> + Debug>(
        path: P,
    ) -> Result<IntoIter<Bytes>, Box<dyn Error>> {
        Ok(Reader::from_file(path)?
            .records()
            .into_iter()
            .par_bridge()
            .map(|read| read.expect("Error reading fasta record."))
            .map(|record| Bytes::copy_from_slice(record.seq()))
            .collect::<Vec<Bytes>>()
            .into_iter())
    }
}