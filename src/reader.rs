use std::{error::Error, fmt::Debug, path::Path};

use bytes::Bytes;
use rayon::{prelude::IntoParallelIterator, vec::IntoIter};

#[cfg(not(feature = "needletail"))]
pub(crate) fn read<P: AsRef<Path> + Debug>(path: P) -> Result<IntoIter<Bytes>, Box<dyn Error>> {
    Ok(bio::io::fasta::Reader::from_file(path)?
        .records()
        .map(|read| read.expect("Error reading FASTA record."))
        .map(|record| Bytes::copy_from_slice(record.seq()))
        .collect::<Vec<Bytes>>()
        .into_par_iter())
}

#[cfg(feature = "needletail")]
pub(crate) fn read<P: AsRef<Path> + Debug>(path: P) -> Result<IntoIter<Bytes>, Box<dyn Error>> {
    let mut reader = needletail::parse_fastx_file(path)?;
    let mut v = Vec::new();
    while let Some(record) = reader.next() {
        let record = record.expect("invalid record");
        let seq = Bytes::copy_from_slice(&record.seq());
        v.push(seq);
    }
    Ok(v.into_par_iter())
}
