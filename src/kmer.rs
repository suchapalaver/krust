use std::cmp::Ordering;

use bytes::Bytes;

#[derive(Debug, Default, Eq, PartialEq, Hash)]
pub struct Kmer {
    pub bytes: Bytes,
    pub packed_bits: u64,
    pub reverse_complement: bool,
    pub count: i32,
}

impl Kmer {
    pub fn from_sub(sub: Bytes) -> Result<Self, usize> {
        sub.into_iter()
            .enumerate()
            .map(|(i, byte)| {
                Ok(match byte {
                    b'A' | b'C' | b'G' | b'T' => byte,
                    _ => return Err(i),
                })
            })
            .collect()
    }

    pub fn pack_bits(&mut self) {
        for elem in self.bytes.iter() {
            self.packed_bits <<= 2;
            let byte: KmerByte = elem.into();
            let mask: u64 = byte.into();
            self.packed_bits |= mask
        }
    }

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

pub enum KmerByte {
    A,
    C,
    G,
    T,
}

impl From<&u8> for KmerByte {
    fn from(val: &u8) -> Self {
        match val {
            b'A' => KmerByte::A,
            b'C' => KmerByte::C,
            b'G' => KmerByte::G,
            b'T' => KmerByte::T,
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
        let sub = &[b'G', b'A', b'T', b'T', b'A', b'C', b'A'];
        let k = Kmer::from_sub(Bytes::copy_from_slice(sub)).unwrap();
        insta::assert_snapshot!(format!("{:?}", k.bytes), @r#"b"GATTACA""#);
    }

    #[test]
    fn from_substring_returns_err_for_invalid_substring() {
        let sub = &[b'N'];
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
}
