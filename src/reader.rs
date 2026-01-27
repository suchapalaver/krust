use std::{error::Error, fmt::Debug, path::Path};

use bytes::Bytes;
use rayon::{prelude::IntoParallelIterator, vec::IntoIter};

use crate::format::SequenceFormat;

/// A sequence with optional quality scores.
///
/// For FASTA input, quality is always `None`.
/// For FASTQ input, quality contains the Phred+33 encoded quality scores.
#[derive(Clone)]
pub struct SequenceWithQuality {
    pub seq: Bytes,
    pub qual: Option<Bytes>,
}

/// Check if a path has a gzip extension (.gz).
#[cfg(all(feature = "gzip", not(feature = "needletail")))]
fn is_gzip_path<P: AsRef<Path>>(path: P) -> bool {
    path.as_ref()
        .extension()
        .map(|ext| ext == "gz")
        .unwrap_or(false)
}

/// Trait for extracting sequence data from records.
#[cfg(not(feature = "needletail"))]
trait SequenceRecord {
    fn seq(&self) -> &[u8];
    fn qual(&self) -> Option<&[u8]>;
}

#[cfg(not(feature = "needletail"))]
impl SequenceRecord for bio::io::fasta::Record {
    fn seq(&self) -> &[u8] {
        bio::io::fasta::Record::seq(self)
    }

    fn qual(&self) -> Option<&[u8]> {
        None
    }
}

#[cfg(not(feature = "needletail"))]
impl SequenceRecord for bio::io::fastq::Record {
    fn seq(&self) -> &[u8] {
        bio::io::fastq::Record::seq(self)
    }

    fn qual(&self) -> Option<&[u8]> {
        Some(bio::io::fastq::Record::qual(self))
    }
}

/// Converts a collection of sequence records to a parallel iterator of Bytes.
#[cfg(not(feature = "needletail"))]
fn records_to_bytes<R: SequenceRecord>(records: Vec<R>) -> IntoIter<Bytes> {
    records
        .into_iter()
        .map(|record| Bytes::copy_from_slice(record.seq()))
        .collect::<Vec<Bytes>>()
        .into_par_iter()
}

/// Converts a collection of sequence records to a parallel iterator with quality.
#[cfg(not(feature = "needletail"))]
fn records_to_bytes_with_quality<R: SequenceRecord>(
    records: Vec<R>,
) -> IntoIter<SequenceWithQuality> {
    records
        .into_iter()
        .map(|record| SequenceWithQuality {
            seq: Bytes::copy_from_slice(record.seq()),
            qual: record.qual().map(Bytes::copy_from_slice),
        })
        .collect::<Vec<SequenceWithQuality>>()
        .into_par_iter()
}

#[cfg(all(not(feature = "needletail"), not(feature = "gzip")))]
pub(crate) fn read<P: AsRef<Path> + Debug>(
    path: P,
    format: SequenceFormat,
) -> Result<IntoIter<Bytes>, Box<dyn Error>> {
    let resolved_format = format.resolve(Some(path.as_ref()));

    match resolved_format {
        SequenceFormat::Fastq => {
            let records: Result<Vec<_>, _> =
                bio::io::fastq::Reader::from_file(path)?.records().collect();
            Ok(records_to_bytes(records?))
        }
        SequenceFormat::Fasta | SequenceFormat::Auto => {
            let records: Result<Vec<_>, _> =
                bio::io::fasta::Reader::from_file(path)?.records().collect();
            Ok(records_to_bytes(records?))
        }
    }
}

#[cfg(all(not(feature = "needletail"), feature = "gzip"))]
pub(crate) fn read<P: AsRef<Path> + Debug>(
    path: P,
    format: SequenceFormat,
) -> Result<IntoIter<Bytes>, Box<dyn Error>> {
    use bio::io::{fasta, fastq};
    use flate2::read::GzDecoder;
    use std::{fs::File, io::BufReader};

    let resolved_format = format.resolve(Some(path.as_ref()));

    if is_gzip_path(&path) {
        let file = File::open(path.as_ref())?;
        let decoder = GzDecoder::new(file);
        let buf_reader = BufReader::new(decoder);

        match resolved_format {
            SequenceFormat::Fastq => {
                let reader = fastq::Reader::new(buf_reader);
                let records: Result<Vec<bio::io::fastq::Record>, _> = reader.records().collect();
                Ok(records_to_bytes(records?))
            }
            SequenceFormat::Fasta | SequenceFormat::Auto => {
                let reader = fasta::Reader::new(buf_reader);
                let records: Result<Vec<bio::io::fasta::Record>, _> = reader.records().collect();
                Ok(records_to_bytes(records?))
            }
        }
    } else {
        match resolved_format {
            SequenceFormat::Fastq => {
                let records: Result<Vec<_>, _> =
                    fastq::Reader::from_file(path)?.records().collect();
                Ok(records_to_bytes(records?))
            }
            SequenceFormat::Fasta | SequenceFormat::Auto => {
                let records: Result<Vec<_>, _> =
                    fasta::Reader::from_file(path)?.records().collect();
                Ok(records_to_bytes(records?))
            }
        }
    }
}

#[cfg(feature = "needletail")]
pub(crate) fn read<P: AsRef<Path> + Debug>(
    path: P,
    _format: SequenceFormat,
) -> Result<IntoIter<Bytes>, Box<dyn Error>> {
    // needletail auto-detects FASTA/FASTQ format via parse_fastx_file
    let mut reader = needletail::parse_fastx_file(path)?;
    let mut v = Vec::new();
    while let Some(record) = reader.next() {
        let record = record?;
        let seq = Bytes::copy_from_slice(&record.seq());
        v.push(seq);
    }
    Ok(v.into_par_iter())
}

/// Reads sequences with quality scores from a file.
///
/// For FASTA input, quality is always `None`.
/// For FASTQ input, quality contains the Phred+33 encoded quality scores.
#[cfg(all(not(feature = "needletail"), not(feature = "gzip")))]
pub(crate) fn read_with_quality<P: AsRef<Path> + Debug>(
    path: P,
    format: SequenceFormat,
) -> Result<IntoIter<SequenceWithQuality>, Box<dyn Error>> {
    let resolved_format = format.resolve(Some(path.as_ref()));

    match resolved_format {
        SequenceFormat::Fastq => {
            let records: Result<Vec<_>, _> =
                bio::io::fastq::Reader::from_file(path)?.records().collect();
            Ok(records_to_bytes_with_quality(records?))
        }
        SequenceFormat::Fasta | SequenceFormat::Auto => {
            let records: Result<Vec<_>, _> =
                bio::io::fasta::Reader::from_file(path)?.records().collect();
            Ok(records_to_bytes_with_quality(records?))
        }
    }
}

/// Reads sequences with quality scores from a file (gzip version).
#[cfg(all(not(feature = "needletail"), feature = "gzip"))]
pub(crate) fn read_with_quality<P: AsRef<Path> + Debug>(
    path: P,
    format: SequenceFormat,
) -> Result<IntoIter<SequenceWithQuality>, Box<dyn Error>> {
    use bio::io::{fasta, fastq};
    use flate2::read::GzDecoder;
    use std::{fs::File, io::BufReader};

    let resolved_format = format.resolve(Some(path.as_ref()));

    if is_gzip_path(&path) {
        let file = File::open(path.as_ref())?;
        let decoder = GzDecoder::new(file);
        let buf_reader = BufReader::new(decoder);

        match resolved_format {
            SequenceFormat::Fastq => {
                let reader = fastq::Reader::new(buf_reader);
                let records: Result<Vec<bio::io::fastq::Record>, _> = reader.records().collect();
                Ok(records_to_bytes_with_quality(records?))
            }
            SequenceFormat::Fasta | SequenceFormat::Auto => {
                let reader = fasta::Reader::new(buf_reader);
                let records: Result<Vec<bio::io::fasta::Record>, _> = reader.records().collect();
                Ok(records_to_bytes_with_quality(records?))
            }
        }
    } else {
        match resolved_format {
            SequenceFormat::Fastq => {
                let records: Result<Vec<_>, _> =
                    fastq::Reader::from_file(path)?.records().collect();
                Ok(records_to_bytes_with_quality(records?))
            }
            SequenceFormat::Fasta | SequenceFormat::Auto => {
                let records: Result<Vec<_>, _> =
                    fasta::Reader::from_file(path)?.records().collect();
                Ok(records_to_bytes_with_quality(records?))
            }
        }
    }
}

/// Reads sequences with quality scores from a file (needletail version).
#[cfg(feature = "needletail")]
pub(crate) fn read_with_quality<P: AsRef<Path> + Debug>(
    path: P,
    _format: SequenceFormat,
) -> Result<IntoIter<SequenceWithQuality>, Box<dyn Error>> {
    let mut reader = needletail::parse_fastx_file(path)?;
    let mut v = Vec::new();
    while let Some(record) = reader.next() {
        let record = record?;
        let seq = Bytes::copy_from_slice(&record.seq());
        let qual = record.qual().map(Bytes::copy_from_slice);
        v.push(SequenceWithQuality { seq, qual });
    }
    Ok(v.into_par_iter())
}
