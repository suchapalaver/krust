use std::cmp::Ordering;

use bytes::Bytes;

#[derive(Debug, Default, Eq, PartialEq, Hash)]
pub(crate) struct Kmer {
    pub(crate) bytes: Bytes,
    pub(crate) packed_bits: u64,
    pub(crate) reverse_complement: bool,
    pub(crate) count: i32,
}

impl Kmer {
    pub(crate) fn from_sub(sub: &Bytes) -> Result<Self, usize> {
        sub.iter()
            .enumerate()
            .map(|x| {
                Ok(match Monomer::try_from(x) {
                    Ok(b) => b,
                    Err(i) => return Err(i),
                })
            })
            .collect()
    }

    pub(crate) fn pack(&mut self) {
        for elem in self.bytes.iter() {
            self.packed_bits <<= 2;
            let mask: u64 = Monomer::from_u8(elem).into();
            self.packed_bits |= mask
        }
    }

    pub(crate) fn canonical(&mut self) {
        match self
            .bytes
            .iter()
            .rev()
            .map(|byte| Monomer::from_u8(byte).complement().into_u8())
            .collect::<Bytes>()
        {
            reverse_complement if reverse_complement.cmp(&self.bytes) == Ordering::Less => {
                self.bytes = reverse_complement;
                self.reverse_complement = true
            }
            _ => (),
        }
    }

    pub(crate) fn unpack(&mut self, k: usize) {
        self.bytes = (0..k)
            .into_iter()
            .map(|i| self.packed_bits << ((i * 2) + 64 - (k * 2)) >> 62)
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

impl TryFrom<(usize, &u8)> for Monomer {
    type Error = usize;

    fn try_from(value: (usize, &u8)) -> Result<Self, Self::Error> {
        match value {
            (_, b'A') => Ok(Self::A),
            (_, b'C') => Ok(Self::C),
            (_, b'G') => Ok(Self::G),
            (_, b'T') => Ok(Self::T),
            (invalid_byte_index, _) => Err(invalid_byte_index),
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

    fn from_u8(byte: &u8) -> Self {
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

#[cfg(test)]
pub mod test {
    use super::*;

    #[test]
    fn bytes_from_valid_substring() {
        let sub = &[b'G', b'A', b'T', b'T', b'A', b'C', b'A'];
        let k = Kmer::from_sub(&Bytes::copy_from_slice(sub)).unwrap();
        insta::assert_snapshot!(format!("{:?}", k.bytes), @r###"b"GATTACA""###);
    }

    #[test]
    fn error_on_invalid_byte() {
        let enumeration = 0;
        let b = b'N';
        assert!(Monomer::try_from((enumeration, &b)).is_err());
    }

    #[test]
    fn from_substring_returns_err_for_invalid_substring() {
        let sub = &[b'N'];
        let k = Kmer::from_sub(&Bytes::copy_from_slice(sub));
        assert!(k.is_err());
    }

    #[test]
    fn from_sub_finds_invalid_byte_index() {
        let dna = "NACNN".as_bytes();
        let ans = Kmer::from_sub(&Bytes::copy_from_slice(dna));
        assert_eq!(Err(0), ans);

        let dna = "ANCNG".as_bytes();
        let ans = Kmer::from_sub(&Bytes::copy_from_slice(dna));
        assert_eq!(Err(1), ans);

        let dna = "AANTG".as_bytes();
        let ans = Kmer::from_sub(&Bytes::copy_from_slice(dna));
        assert_eq!(Err(2), ans);

        let dna = "CCCNG".as_bytes();
        let ans = Kmer::from_sub(&Bytes::copy_from_slice(dna));
        assert_eq!(Err(3), ans);

        let dna = "AACTN".as_bytes();
        let ans = Kmer::from_sub(&Bytes::copy_from_slice(dna));
        assert_eq!(Err(4), ans);
    }
}
