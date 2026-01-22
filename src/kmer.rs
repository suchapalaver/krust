//! K-mer representation and manipulation.
//!
//! This module provides types for working with DNA k-mers, including:
//! - Bit-packing k-mers into 64-bit integers (supports k <= 32)
//! - Computing canonical k-mers (lexicographically smaller of k-mer/reverse-complement)
//! - Converting between byte and packed representations

use std::cmp::Ordering;

use bytes::Bytes;

/// A DNA k-mer with both byte and packed bit representations.
///
/// K-mers are short DNA sequences of length k. This struct stores both the
/// raw byte representation and a packed 64-bit integer representation for
/// efficient hashing and comparison.
#[derive(Debug, Default, Clone, Eq, PartialEq, Hash)]
pub struct Kmer {
    /// The k-mer as ASCII bytes (A, C, G, T).
    pub bytes: Bytes,
    /// The k-mer packed as a 64-bit integer (2 bits per base).
    pub packed_bits: u64,
    /// Whether this k-mer is the reverse complement of the original.
    pub reverse_complement: bool,
    /// The count/frequency of this k-mer.
    pub count: i32,
}

impl Kmer {
    /// Creates a k-mer from a byte sequence.
    ///
    /// Returns `Ok(Kmer)` if all bytes are valid DNA bases (A, C, G, T, or lowercase a, c, g, t).
    /// Soft-masked (lowercase) bases are converted to uppercase.
    /// Returns `Err(index)` with the index of the first invalid byte.
    pub fn from_sub(sub: Bytes) -> Result<Self, usize> {
        sub.into_iter()
            .enumerate()
            .map(|(i, byte)| {
                Ok(match byte {
                    b'A' | b'C' | b'G' | b'T' => byte,
                    b'a' | b'c' | b'g' | b't' => byte.to_ascii_uppercase(),
                    _ => return Err(i),
                })
            })
            .collect()
    }

    /// Packs the k-mer bytes into a 64-bit integer.
    ///
    /// Each base is encoded as 2 bits: A=00, C=01, G=10, T=11.
    /// This allows k-mers up to length 32 to fit in a u64.
    pub fn pack_bits(&mut self) {
        for elem in self.bytes.iter() {
            self.packed_bits <<= 2;
            let byte: KmerByte = elem.into();
            let mask: u64 = byte.into();
            self.packed_bits |= mask
        }
    }

    /// Converts the k-mer to its canonical form.
    ///
    /// The canonical form is the lexicographically smaller of the k-mer
    /// and its reverse complement. Sets `reverse_complement` to true if
    /// the reverse complement was chosen.
    pub fn canonical(&mut self) {
        match self
            .bytes
            .iter()
            .rev()
            .map(KmerByte::from)
            .map(KmerByte::reverse_complement)
            .map(KmerByte::into)
            .collect::<Bytes>()
        {
            reverse_complement if reverse_complement.cmp(&self.bytes) == Ordering::Less => {
                self.bytes = reverse_complement;
                self.reverse_complement = true
            }
            _ => (),
        }
    }

    /// Unpacks a 64-bit integer back into k-mer bytes.
    ///
    /// The `k` parameter specifies the k-mer length to extract.
    pub fn unpack_bits(&mut self, k: usize) {
        self.bytes = (0..k)
            .map(|i| self.packed_bits << ((i * 2) + 64 - (k * 2)) >> 62)
            .map(KmerByte::from)
            .map(KmerByte::into)
            .collect()
    }
}

impl FromIterator<u8> for Kmer {
    fn from_iter<I: IntoIterator<Item = u8>>(iter: I) -> Self {
        Self {
            bytes: iter.into_iter().collect(),
            ..Default::default()
        }
    }
}

/// A single DNA base (nucleotide).
///
/// Used for converting between byte representations (ASCII) and
/// numeric representations (for bit packing).
pub enum KmerByte {
    /// Adenine
    A,
    /// Cytosine
    C,
    /// Guanine
    G,
    /// Thymine
    T,
}

impl From<&u8> for KmerByte {
    fn from(val: &u8) -> Self {
        match val {
            b'A' | b'a' => KmerByte::A,
            b'C' | b'c' => KmerByte::C,
            b'G' | b'g' => KmerByte::G,
            b'T' | b't' => KmerByte::T,
            _ => unreachable!(),
        }
    }
}

impl From<KmerByte> for u8 {
    fn from(val: KmerByte) -> Self {
        match val {
            KmerByte::A => b'A',
            KmerByte::C => b'C',
            KmerByte::G => b'G',
            KmerByte::T => b'T',
        }
    }
}

impl From<u64> for KmerByte {
    fn from(val: u64) -> Self {
        match val {
            0 => KmerByte::A,
            1 => KmerByte::C,
            2 => KmerByte::G,
            3 => KmerByte::T,
            _ => unreachable!(),
        }
    }
}

impl From<KmerByte> for u64 {
    fn from(val: KmerByte) -> Self {
        match val {
            KmerByte::A => 0,
            KmerByte::C => 1,
            KmerByte::G => 2,
            KmerByte::T => 3,
        }
    }
}

impl KmerByte {
    pub fn reverse_complement(self) -> Self {
        match self {
            Self::A => Self::T,
            Self::C => Self::G,
            Self::G => Self::C,
            Self::T => Self::A,
        }
    }
}

#[cfg(test)]
pub mod test {
    use super::*;

    #[test]
    fn bytes_from_valid_substring() {
        let sub = b"GATTACA";
        let k = Kmer::from_sub(Bytes::copy_from_slice(sub)).unwrap();
        insta::assert_snapshot!(format!("{:?}", k.bytes), @r#"b"GATTACA""#);
    }

    #[test]
    fn from_substring_returns_err_for_invalid_substring() {
        let sub = b"N";
        let k = Kmer::from_sub(Bytes::copy_from_slice(sub));
        assert!(k.is_err());
    }

    #[test]
    fn from_sub_finds_invalid_byte_index() {
        let dna = "NACNN".as_bytes();
        let res = Kmer::from_sub(Bytes::copy_from_slice(dna));
        assert_eq!(Err(0), res);

        let dna = "ANCNG".as_bytes();
        let res = Kmer::from_sub(Bytes::copy_from_slice(dna));
        assert_eq!(Err(1), res);

        let dna = "AANTG".as_bytes();
        let res = Kmer::from_sub(Bytes::copy_from_slice(dna));
        assert_eq!(Err(2), res);

        let dna = "CCCNG".as_bytes();
        let res = Kmer::from_sub(Bytes::copy_from_slice(dna));
        assert_eq!(Err(3), res);

        let dna = "AACTN".as_bytes();
        let res = Kmer::from_sub(Bytes::copy_from_slice(dna));
        assert_eq!(Err(4), res);
    }

    #[test]
    fn pack_unpack_roundtrip() {
        let sequences = ["ACGT", "AAAA", "TTTT", "CCCC", "GGGG", "GATTACA"];
        for seq in sequences {
            let mut kmer = Kmer::from_sub(Bytes::copy_from_slice(seq.as_bytes())).unwrap();
            kmer.pack_bits();
            kmer.unpack_bits(seq.len());
            assert_eq!(kmer.bytes.as_ref(), seq.as_bytes());
        }
    }

    #[test]
    fn pack_unpack_roundtrip_various_lengths() {
        for k in 1..=32 {
            let seq = "A".repeat(k);
            let mut kmer = Kmer::from_sub(Bytes::copy_from_slice(seq.as_bytes())).unwrap();
            kmer.pack_bits();
            kmer.unpack_bits(k);
            assert_eq!(kmer.bytes.as_ref(), seq.as_bytes());
        }
    }

    #[test]
    fn canonical_selects_lexicographically_smaller() {
        // ACGT and ACGT (reverse complement) - ACGT is palindromic
        let mut kmer = Kmer::from_sub(Bytes::copy_from_slice(b"ACGT")).unwrap();
        kmer.canonical();
        assert_eq!(kmer.bytes.as_ref(), b"ACGT");
        assert!(!kmer.reverse_complement);

        // AAA -> TTT reverse complement, AAA is smaller
        let mut kmer = Kmer::from_sub(Bytes::copy_from_slice(b"AAA")).unwrap();
        kmer.canonical();
        assert_eq!(kmer.bytes.as_ref(), b"AAA");
        assert!(!kmer.reverse_complement);

        // TTT -> AAA reverse complement, AAA is smaller
        let mut kmer = Kmer::from_sub(Bytes::copy_from_slice(b"TTT")).unwrap();
        kmer.canonical();
        assert_eq!(kmer.bytes.as_ref(), b"AAA");
        assert!(kmer.reverse_complement);

        // GATTACA -> TGTAATC reverse complement, GATTACA is smaller
        let mut kmer = Kmer::from_sub(Bytes::copy_from_slice(b"GATTACA")).unwrap();
        kmer.canonical();
        assert_eq!(kmer.bytes.as_ref(), b"GATTACA");
        assert!(!kmer.reverse_complement);

        // TGTAATC -> GATTACA reverse complement, GATTACA is smaller
        let mut kmer = Kmer::from_sub(Bytes::copy_from_slice(b"TGTAATC")).unwrap();
        kmer.canonical();
        assert_eq!(kmer.bytes.as_ref(), b"GATTACA");
        assert!(kmer.reverse_complement);
    }

    #[test]
    fn kmer_byte_reverse_complement() {
        assert!(matches!(KmerByte::A.reverse_complement(), KmerByte::T));
        assert!(matches!(KmerByte::T.reverse_complement(), KmerByte::A));
        assert!(matches!(KmerByte::C.reverse_complement(), KmerByte::G));
        assert!(matches!(KmerByte::G.reverse_complement(), KmerByte::C));
    }

    #[test]
    fn kmer_byte_to_u64_roundtrip() {
        for (byte, expected) in [(b'A', 0u64), (b'C', 1u64), (b'G', 2u64), (b'T', 3u64)] {
            let kmer_byte: KmerByte = (&byte).into();
            let val: u64 = kmer_byte.into();
            assert_eq!(val, expected);
            let back: KmerByte = val.into();
            let back_byte: u8 = back.into();
            assert_eq!(back_byte, byte);
        }
    }

    #[test]
    fn empty_kmer() {
        let kmer = Kmer::from_sub(Bytes::new()).unwrap();
        assert!(kmer.bytes.is_empty());
    }

    #[test]
    fn single_base_kmers() {
        for base in [b'A', b'C', b'G', b'T'] {
            let kmer = Kmer::from_sub(Bytes::copy_from_slice(&[base])).unwrap();
            assert_eq!(kmer.bytes.as_ref(), &[base]);
        }
    }

    #[test]
    fn soft_masked_bases_converted_to_uppercase() {
        // Test that lowercase bases (soft-masked) are accepted and converted
        let sub = b"AAAa";
        let k = Kmer::from_sub(Bytes::copy_from_slice(sub)).unwrap();
        insta::assert_snapshot!(format!("{:?}", k.bytes), @r#"b"AAAA""#);
    }

    #[test]
    fn soft_masked_all_lowercase() {
        // Test all lowercase sequence
        let sub = b"gattaca";
        let k = Kmer::from_sub(Bytes::copy_from_slice(sub)).unwrap();
        insta::assert_snapshot!(format!("{:?}", k.bytes), @r#"b"GATTACA""#);
    }

    #[test]
    fn soft_masked_mixed_case() {
        // Test mixed case sequence
        let sub = b"AcGt";
        let k = Kmer::from_sub(Bytes::copy_from_slice(sub)).unwrap();
        insta::assert_snapshot!(format!("{:?}", k.bytes), @r#"b"ACGT""#);
    }
}
