use bytes::Bytes;

custom_error::custom_error! { pub ValidityError
    InvalidByte = "not a valid byte",
}

/// A valid k-mer
#[derive(Debug, Eq, PartialEq)]
pub struct Kmer(pub Bytes);

impl Kmer {
    pub(crate) fn from_sub(sub: &Bytes) -> Result<Self, ValidityError> {
        sub.iter().map(|b| Monomer::try_from(*b)).collect()
    }

    pub(crate) fn canonical<'a>(
        reverse_complement: &'a Bytes,
        kmer: &'a Bytes,
    ) -> &'a Bytes {
        match reverse_complement.cmp(kmer) {
            std::cmp::Ordering::Less => reverse_complement,
            _ => kmer,
        }
    }

    pub(crate) fn find_invalid(sub: &Bytes) -> usize {
        sub.iter()
            .rposition(|b| Monomer::try_from(*b).is_err())
            .unwrap()
    }

    pub(crate) fn bytes(&self) -> &Bytes {
        &self.0
    }
}

impl FromIterator<u8> for Kmer {
    fn from_iter<I: IntoIterator<Item = u8>>(iter: I) -> Self {
        Self(iter.into_iter().collect())
    }
}

pub(crate) enum Monomer {
    A,
    C,
    G,
    T,
}

impl TryFrom<u8> for Monomer {
    type Error = ValidityError;

    fn try_from(value: u8) -> Result<Self, Self::Error> {
        match value {
            b'A' => Ok(Self::A),
            b'C' => Ok(Self::C),
            b'G' => Ok(Self::G),
            b'T' => Ok(Self::T),
            _ => Err(ValidityError::InvalidByte),
        }
    }
}

impl FromIterator<Monomer> for Kmer {
    fn from_iter<I: IntoIterator<Item = Monomer>>(iter: I) -> Self {
        Self(iter.into_iter().map(Monomer::into_u8).collect())
    }
}

impl From<u64> for Monomer {
    fn from(u: u64) -> Self {
        match u {
            0 => Self::A,
            1 => Self::C,
            2 => Self::G,
            _ => Self::T,
        }
    }
}

#[allow(clippy::from_over_into)]
impl Into<u64> for Monomer {
    fn into(self) -> u64 {
        match self {
            Self::A => 0,
            Self::C => 1,
            Self::G => 2,
            Self::T => 3,
        }
    }
}

impl Monomer {
    fn complement(self) -> Self {
        match self {
            Self::A => Self::T,
            Self::C => Self::G,
            Self::G => Self::C,
            Self::T => Self::A,
        }
    }

    fn from_u8(byte: u8) -> Self {
        match byte {
            b'A' => Self::A,
            b'C' => Self::C,
            b'G' => Self::G,
            _ => Self::T,
        }
    }

    fn into_u8(self) -> u8 {
        match self {
            Self::A => b'A',
            Self::C => b'C',
            Self::G => b'G',
            Self::T => b'T',
        }
    }
}

trait Pack {
    fn isolate(self, i: usize, k: usize) -> Self;
    fn replace(self) -> Self;
}

impl Pack for u64 {
    fn isolate(self, i: usize, k: usize) -> Self {
        self << ((i * 2) + 64 - (k * 2))
    }

    fn replace(self) -> Self {
        self >> 62
    }
}

/// Compressing k-mers of length `0 < k < 33`, bitpacking them into unsigned integers
pub(crate) struct Bitpack(pub u64);

impl Bitpack {
    fn new() -> Self {
        Self(0)
    }

    fn pack(&mut self, elem: &u8) {
        self.shift();
        let mask: u64 = Monomer::from_u8(*elem).into();
        self.0 |= mask
    }

    fn shift(&mut self) {
        self.0 <<= 2
    }
}

impl From<&Bytes> for Bitpack {
    fn from(bytes: &Bytes) -> Self {
        let mut packed = Self::new();
        bytes.iter().for_each(|b| {
            packed.pack(b)
        });
        packed
    }
}

impl From<&Kmer> for Bitpack {
    fn from(kmer: &Kmer) -> Self {
        let mut packed = Self::new();
        kmer.bytes().iter().for_each(|b| {
            packed.pack(b)
        });
        packed
    }
}

/// Unpack bitpacked k-mer data
#[derive(Hash, PartialEq, Eq)]
pub struct Unpack(pub Bytes);

impl Unpack {
    pub(crate) fn bit(bit: u64, k: usize) -> Self {
        (0..k)
            .into_iter()
            .map(|i| bit.isolate(i, k))
            .map(|bit| bit.replace())
            .map(Monomer::from)
            .map(|m| m.into_u8())
            .collect()
    }
}

impl FromIterator<u8> for Unpack {
    fn from_iter<I: IntoIterator<Item = u8>>(iter: I) -> Self {
        Self(iter.into_iter().collect())
    }
}

/// Convert a DNA string slice into its [reverse compliment](https://en.wikipedia.org/wiki/Complementarity_(molecular_biology)#DNA_and_RNA_base_pair_complementarity).
pub(crate) struct RevComp(pub Bytes);

impl FromIterator<u8> for RevComp {
    fn from_iter<I: IntoIterator<Item = u8>>(iter: I) -> Self {
        Self(iter.into_iter().collect())
    }
}

impl RevComp {
    pub(crate) fn from_kmer(kmer: &Kmer) -> Self {
        kmer.0
            .iter()
            .rev()
            .map(|byte| Monomer::from_u8(*byte).complement().into_u8())
            .collect()
    }
}

#[cfg(test)]
pub mod test {
    use super::*;

    #[test]
    fn test_from_valid_substring() {
        let sub = &[b'G', b'A', b'T', b'T', b'A', b'C', b'A'];
        let k = Kmer::from_sub(&Bytes::copy_from_slice(sub)).unwrap();
        insta::assert_snapshot!(format!("{:?}", k), @r###"Kmer(b"GATTACA")"###);
    }

    #[test]
    fn test_parse_valid_byte() {
        let b = b'N';
        assert!(Monomer::try_from(b).is_err());
    }

    #[test]
    fn from_substring_returns_err_for_invalid_substring() {
        let sub = &[b'N'];
        let k = Kmer::from_sub(&Bytes::copy_from_slice(sub));
        assert!(k.is_err());
    }

    #[test]
    fn find_invalid_works() {
        let dna = "NACNN".as_bytes();
        let ans = Kmer::find_invalid(&Bytes::copy_from_slice(dna));
        assert_eq!(4, ans);
        assert_eq!(&b'N', dna.iter().collect::<Vec<_>>()[ans]);

        let dna = "NACNG".as_bytes();
        let ans = Kmer::find_invalid(&Bytes::copy_from_slice(dna));
        assert_eq!(3, ans);
        assert_eq!(&b'N', dna.iter().collect::<Vec<_>>()[ans]);

        let dna = "NANTG".as_bytes();
        let ans = Kmer::find_invalid(&Bytes::copy_from_slice(dna));
        assert_eq!(2, ans);
        assert_eq!(&b'N', dna.iter().collect::<Vec<_>>()[ans]);

        let dna = "NNCTG".as_bytes();
        let ans = Kmer::find_invalid(&Bytes::copy_from_slice(dna));
        assert_eq!(1, ans);
        assert_eq!(&b'N', dna.iter().collect::<Vec<_>>()[ans]);

        let dna = "NACTG".as_bytes();
        let ans = Kmer::find_invalid(&Bytes::copy_from_slice(dna));
        assert_eq!(0, ans);
        assert_eq!(&b'N', dna.iter().collect::<Vec<_>>()[ans]);
    }
}
