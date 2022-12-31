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
    pub(crate) fn from_sub(sub: Bytes) -> Result<Self, usize> {
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

    pub(crate) fn pack(&mut self) {
        for elem in self.bytes.iter() {
            self.packed_bits <<= 2;
            let mask: u64 = u8::into_u64(elem);
            self.packed_bits |= mask
        }
    }

    pub(crate) fn canonical(&mut self) {
        match self
            .bytes
            .iter()
            .rev()
            .map(u8::complement)
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
            .map(u8::from_u64)
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

trait ByteConversions {
    fn complement(u: &Self) -> Self;
    fn from_u64(u: u64) -> Self;
    fn into_u64(u: &Self) -> u64;
}

impl ByteConversions for u8 {
    fn complement(u: &Self) -> Self {
        match *u {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            _ => b'A',
        }
    }

    fn from_u64(u: u64) -> Self {
        match u {
            0 => b'A',
            1 => b'C',
            2 => b'G',
            _ => b'T',
        }
    }

    fn into_u64(u: &Self) -> u64 {
        match *u {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            _ => 3,
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
        insta::assert_snapshot!(format!("{:?}", k.bytes), @r###"b"GATTACA""###);
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
        let ans = Kmer::from_sub(Bytes::copy_from_slice(dna));
        assert_eq!(Err(0), ans);

        let dna = "ANCNG".as_bytes();
        let ans = Kmer::from_sub(Bytes::copy_from_slice(dna));
        assert_eq!(Err(1), ans);

        let dna = "AANTG".as_bytes();
        let ans = Kmer::from_sub(Bytes::copy_from_slice(dna));
        assert_eq!(Err(2), ans);

        let dna = "CCCNG".as_bytes();
        let ans = Kmer::from_sub(Bytes::copy_from_slice(dna));
        assert_eq!(Err(3), ans);

        let dna = "AACTN".as_bytes();
        let ans = Kmer::from_sub(Bytes::copy_from_slice(dna));
        assert_eq!(Err(4), ans);
    }
}
