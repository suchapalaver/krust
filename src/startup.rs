use bio::io::fasta;
use dashmap::DashMap;
use fxhash::FxHasher;
use rayon::prelude::*;
use std::{
    collections::HashMap,
    error::Error,
    hash::BuildHasherDefault,
    io::{BufWriter, Stdout, Write},
};

pub fn run(filepath: String, k: usize) -> Result<(), Box<dyn Error>> {
    let mut buf = BufWriter::new(std::io::stdout());

    let _print_results = build_kmer_map(filepath, k)?
        .into_iter()
        .par_bridge()
        .map(|(bitpacked_kmer, freq)| (UnpackedKmer::from((bitpacked_kmer, k)).0, freq))
        .map(|(unpacked_kmer, freq)| {
            let kmer_str = String::from_utf8(unpacked_kmer).unwrap();
            (kmer_str, freq)
        })
        .collect::<HashMap<String, i32>>()
        .into_iter()
        .for_each(|(kmer, count)| {
            print_kmer_map(&mut buf, kmer, count);
        });

    buf.flush()?;

    Ok(())
}

/// A custom `DashMap` w/ `FxHasher`.
///
/// ```use dashmap::DashMap;```
/// ```use fxhash::FxHasher;```
/// ```// skip```
/// ```let dashfx_hash: DashFx = DashMap::with_hasher(BuildHasherDefault::<FxHasher>::default());```
/// Useful: [Using a Custom Hash Function in Rust](https://docs.rs/hashers/1.0.1/hashers/#using-a-custom-hash-function-in-rust).
pub type DashFx = DashMap<u64, i32, BuildHasherDefault<FxHasher>>;

pub type UnpackDashFx = DashMap<Vec<u8>, i32, BuildHasherDefault<FxHasher>>;

///  - Reads sequences from fasta records in parallel using [`rayon`](https://docs.rs/rayon/1.5.1/rayon/),
/// using a customized [`dashmap`](https://docs.rs/dashmap/4.0.2/dashmap/struct.DashMap.html)
/// with [`FxHasher`](https://docs.rs/fxhash/0.2.1/fxhash/struct.FxHasher.html) to update in parallel a
/// hashmap of canonical k-mers (keys) and their frequency in the data (values).
fn build_kmer_map(filepath: String, k: usize) -> Result<DashFx, Box<dyn Error>> {
    let kmer_map: DashFx = DashMap::with_hasher(BuildHasherDefault::<FxHasher>::default());
    let _ = fasta::Reader::from_file(&filepath)?
        .records()
        .into_iter()
        .par_bridge()
        .for_each(|r| {
            let record = r.expect("Error reading fasta record.");
            let seq: &[u8] = record.seq();
            process_seq(seq, &k, &kmer_map);
        });
    Ok(kmer_map)
}

///  - Ignores substrings containing `N`.
///  - Canonicalizes by lexicographically smaller of k-mer/reverse-complement.
fn process_seq(seq: &[u8], k: &usize, kmer_map: &DashFx) {
    let mut i = 0;
    while i <= seq.len() - k {
        let sub = &seq[i..i + k];
        let bytestring = Kmer::new(sub);
        if let Some(Kmer(valid_bytestring)) = bytestring {
            process_valid_bytes(kmer_map, valid_bytestring);
            i += 1;
        } else {
            let invalid_byte_index = Kmer::find_invalid(sub);
            i += invalid_byte_index + 1;
        }
    }
}

/// Converts a valid sequence substring from a bytes string to a u64.
fn process_valid_bytes(kmer_map: &DashFx, valid_bytestring: Vec<u8>) {
    let BitpackedKmer(bitpacked_kmer) = BitpackedKmer::from(&valid_bytestring);
    // If the k-mer as found in the sequence is already a key in the `Dashmap`,
    // increment its value and move on.
    if let Some(mut freq) = kmer_map.get_mut(&bitpacked_kmer) {
        *freq += 1;
    } else {
        // Initialize the reverse complement of this so-far unrecorded k-mer.
        let RevCompKmer(revcompkmer) = RevCompKmer::from(&valid_bytestring);
        // Find the alphabetically less of the k-mer substring and its reverse complement.
        let CanonicalKmer(canonical_kmer) = CanonicalKmer::from((revcompkmer, valid_bytestring));
        // Compress the canonical k-mer into a bitpacked 64-bit unsigned integer.
        let BitpackedKmer(kmer) = BitpackedKmer::from(&canonical_kmer);
        // Add k-mer key and initial value to results.
        *kmer_map.entry(kmer).or_insert(0) += 1;
    }
}

fn print_kmer_map(buf: &mut BufWriter<Stdout>, kmer: String, count: i32) {
    writeln!(buf, ">{}\n{}", count, kmer).expect("Unable to write output.");
}

/// Creating a valid k-mer bytestring.
#[derive(Debug, PartialEq)]
pub struct Kmer(Vec<u8>);

impl Kmer {
    pub fn new(sub: &[u8]) -> Option<Kmer> {
        if !sub.contains(&b'N') {
            let valid_kmer = sub.to_vec();
            Some(Kmer(valid_kmer))
        } else {
            None
        }
    }

    /// Find the index of the rightmost invalid byte in an invalid bytestring.
    pub fn find_invalid(sub: &[u8]) -> usize {
        match sub
            .iter()
            .rposition(|byte| ![b'A', b'C', b'G', b'T'].contains(byte))
        {
            Some(rightmost_invalid_byte_index) => rightmost_invalid_byte_index,
            None => panic!("Valid bytestring passed to `find_invalid`, which is a bug."),
        }
    }
}

/// Compressing k-mers of length `0 < k < 33`, bitpacking them into unsigned integers.
pub struct BitpackedKmer(u64);

impl From<&Vec<u8>> for BitpackedKmer {
    fn from(sub: &Vec<u8>) -> Self {
        let bitpacked_kmer: u64 = {
            let mut k: u64 = 0;
            for byte in sub.iter() {
                k <<= 2;
                let mask = match *byte {
                    b'A' => 0,
                    b'C' => 1,
                    b'G' => 2,
                    b'T' => 3,
                    _ => panic!("`BitpackerKmer` handling an invalid k-mer bytestring is unexpected behavior"),
                };
                k |= mask;
            }
            k
        };
        BitpackedKmer(bitpacked_kmer)
    }
}

/// Converting a DNA string slice into its [reverse compliment](https://en.wikipedia.org/wiki/Complementarity_(molecular_biology)#DNA_and_RNA_base_pair_complementarity).
pub struct RevCompKmer(Vec<u8>);

impl From<&Vec<u8>> for RevCompKmer {
    fn from(sub: &Vec<u8>) -> Self {
        let mut revcomp = Vec::with_capacity(sub.len());

        for byte in sub.iter().rev() {
            let comp = RevCompKmer::complement(*byte);
            revcomp.push(comp);
        }
        RevCompKmer(revcomp)
    }
}

impl RevCompKmer {
    fn complement(byte: u8) -> u8 {
        match byte {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            _ => panic!("`RevCompKmer::from` should only be passed valid k-mers"),
        }
    }
}

/// Find the canonical kmer
/// --the alphabetically smaller of the substring and its reverse complement.
pub struct CanonicalKmer(Vec<u8>);

impl From<(Vec<u8>, Vec<u8>)> for CanonicalKmer {
    fn from(comp: (Vec<u8>, Vec<u8>)) -> Self {
        let (reverse_complement, kmer) = (comp.0, comp.1);
        let canonical_kmer = if reverse_complement.cmp(&kmer) == std::cmp::Ordering::Less {
            reverse_complement
        } else {
            kmer
        };
        CanonicalKmer(canonical_kmer)
    }
}

/// Unpacking compressed, bitpacked k-mer data.
#[derive(Hash, PartialEq, Eq)]
pub struct UnpackedKmer(Vec<u8>);

impl From<(u64, usize)> for UnpackedKmer {
    fn from(kmer_data: (u64, usize)) -> Self {
        let (kmer, k) = (kmer_data.0, kmer_data.1);
        let mut byte_string = Vec::with_capacity(k);
        for i in 0..k {
            let isolate = kmer << ((i * 2) + 64 - (k * 2));
            let base = isolate >> 62;
            let byte = UnpackedKmerByte::from(base);
            byte_string.push(byte.0);
        }
        UnpackedKmer(byte_string)
    }
}

/// Unpacking compressed, bitpacked k-mer data.
pub struct UnpackedKmerByte(u8);

impl From<u64> for UnpackedKmerByte {
    fn from(base: u64) -> Self {
        let unpacked_byte = match base {
            0 => b'A',
            1 => b'C',
            2 => b'G',
            3 => b'T',
            _ => panic!("An invalid k-mer passed to here means we have a serious bug"),
        };
        UnpackedKmerByte(unpacked_byte)
    }
}
