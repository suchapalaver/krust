use std::{error::Error, fmt::Debug, path::Path};

use bytes::Bytes;
use rayon::{prelude::IntoParallelIterator, vec::IntoIter};

/// Check if a path has a gzip extension (.gz).
#[cfg(all(feature = "gzip", not(feature = "needletail")))]
fn is_gzip_path<P: AsRef<Path>>(path: P) -> bool {
    path.as_ref()
        .extension()
        .map(|ext| ext == "gz")
        .unwrap_or(false)
}

#[cfg(all(not(feature = "needletail"), not(feature = "gzip")))]
pub(crate) fn read<P: AsRef<Path> + Debug>(path: P) -> Result<IntoIter<Bytes>, Box<dyn Error>> {
    let records: Result<Vec<_>, _> = bio::io::fasta::Reader::from_file(path)?.records().collect();
    Ok(records?
        .into_iter()
        .map(|record| Bytes::copy_from_slice(record.seq()))
        .collect::<Vec<Bytes>>()
        .into_par_iter())
}

#[cfg(all(not(feature = "needletail"), feature = "gzip"))]
pub(crate) fn read<P: AsRef<Path> + Debug>(path: P) -> Result<IntoIter<Bytes>, Box<dyn Error>> {
    use bio::io::fasta;
    use flate2::read::GzDecoder;
    use std::{fs::File, io::BufReader};

    if is_gzip_path(&path) {
        let file = File::open(path.as_ref())?;
        let decoder = GzDecoder::new(file);
        let reader = fasta::Reader::new(BufReader::new(decoder));
        let records: Result<Vec<bio::io::fasta::Record>, _> = reader.records().collect();
        Ok(records?
            .into_iter()
            .map(|record| Bytes::copy_from_slice(record.seq()))
            .collect::<Vec<Bytes>>()
            .into_par_iter())
    } else {
        let records: Result<Vec<_>, _> = fasta::Reader::from_file(path)?.records().collect();
        Ok(records?
            .into_iter()
            .map(|record| Bytes::copy_from_slice(record.seq()))
            .collect::<Vec<Bytes>>()
            .into_par_iter())
    }
}

#[cfg(feature = "needletail")]
pub(crate) fn read<P: AsRef<Path> + Debug>(path: P) -> Result<IntoIter<Bytes>, Box<dyn Error>> {
    let mut reader = needletail::parse_fastx_file(path)?;
    let mut v = Vec::new();
    while let Some(record) = reader.next() {
        let record = record?;
        let seq = Bytes::copy_from_slice(&record.seq());
        v.push(seq);
    }
    Ok(v.into_par_iter())
}
