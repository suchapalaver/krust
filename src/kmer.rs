use std::cmp::Ordering;

use bytes::Bytes;

custom_error::custom_error! { pub ValidityError
    InvalidByte = "not a valid byte",
}

/// Compressing k-mers of length `0 < k < 33`, Kmering them into unsigned integers
#[derive(Debug, Default, Eq, PartialEq, Hash)]
pub(crate) struct Kmer {
    pub(crate) bytes: Bytes,
    pub(crate) reverse_complement: Bytes,
    pub(crate) packed_bits: u64,
    pub(crate) count: i32,
}

impl Kmer {
    pub(crate) fn from_sub(sub: &Bytes) -> Result<Self, ValidityError> {
        sub.iter().map(|b| Monomer::try_from(*b)).collect()
    }

    pub(crate) fn find_invalid(sub: &Bytes) -> usize {
        sub.iter()
            .rposition(|b| Monomer::try_from(*b).is_err())
            .unwrap()
    }

    pub(crate) fn pack(&mut self) {
        let iter = if self.bytes.is_empty() {
            &self.reverse_complement
        } else {
            &self.bytes
        };
        for elem in iter.iter() {
            self.packed_bits <<= 2;
            let mask: u64 = Monomer::from_u8(*elem).into();
            self.packed_bits |= mask
        }
    }

    pub(crate) fn reverse_complement(&mut self) {
        self.reverse_complement = self
            .bytes
            .iter()
            .rev()
            .map(|byte| Monomer::from_u8(*byte).complement().into_u8())
            .collect();
    }

    pub(crate) fn canonical(&mut self) {
        match self.reverse_complement.cmp(&self.bytes) {
            Ordering::Less => self.bytes.clear(),
            _ => self.reverse_complement.clear(),
        };
    }

    pub(crate) fn unpack(&mut self, k: usize) {
        self.bytes = (0..k)
            .into_iter()
            .map(|i| self.packed_bits.isolate(i, k))
            .map(|bit| bit.replace())
            .map(Monomer::from)
            .map(|m| m.into_u8())
            .collect()
    }
}

impl FromIterator<Monomer> for Kmer {
    fn from_iter<I: IntoIterator<Item = Monomer>>(iter: I) -> Self {
        Self {
            bytes: iter.into_iter().map(Monomer::into_u8).collect(),
            ..Default::default()
        }
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

impl From<Monomer> for u64 {
    fn from(m: Monomer) -> u64 {
        match m {
            Monomer::A => 0,
            Monomer::C => 1,
            Monomer::G => 2,
            Monomer::T => 3,
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

#[cfg(test)]
pub mod test {
    use super::*;

    #[test]
    fn test_from_valid_substring() {
        let sub = &[b'G', b'A', b'T', b'T', b'A', b'C', b'A'];
        let k = Kmer::from_sub(&Bytes::copy_from_slice(sub)).unwrap();
        insta::assert_snapshot!(format!("{:?}", k.bytes), @r###"b"GATTACA""###);
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
