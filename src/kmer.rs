custom_error::custom_error! { pub ValidityError
    InvalidByte = "not a valid byte",
}

/// A valid k-mer
#[derive(Debug, Eq, PartialEq)]
pub struct Kmer(pub Vec<u8>);

impl Kmer {
    fn new() -> Self {
        Self(Vec::new())
    }

    fn add(&mut self, elem: u8) {
        self.0.push(elem)
    }

    pub(crate) fn from_sub(sub: &[u8]) -> Result<Self, ValidityError> {
        sub.iter().cloned().map(Monomer::parse).collect()
    }

    pub(crate) fn canonical(
        reverse_complement: Vec<u8>,
        kmer: Vec<u8>,
    ) -> impl Iterator<Item = u8> {
        match reverse_complement.cmp(&kmer) {
            std::cmp::Ordering::Less => reverse_complement.into_iter(),
            _ => kmer.into_iter(),
        }
    }

    pub(crate) fn find_invalid(sub: &[u8]) -> usize {
        sub.iter()
            .rposition(|byte| Monomer::parse(*byte).is_err())
            .unwrap()
    }
}

impl FromIterator<u8> for Kmer {
    fn from_iter<I: IntoIterator<Item = u8>>(iter: I) -> Self {
        let mut k = Kmer::new();

        for i in iter {
            k.add(i);
        }
        k
    }
}

pub(crate) enum Monomer {
    A(u8),
    C(u8),
    G(u8),
    T(u8),
}

impl From<u8> for Monomer {
    fn from(m: u8) -> Self {
        match m {
            b'A' => Self::A(m),
            b'C' => Self::C(m),
            b'G' => Self::G(m),
            _ => Self::T(m),
        }
    }
}

impl From<Monomer> for u8 {
    fn from(m: Monomer) -> Self {
        match m {
            Monomer::A(b) => b,
            Monomer::C(b) => b,
            Monomer::G(b) => b,
            Monomer::T(b) => b,
        }
    }
}

impl FromIterator<Monomer> for Kmer {
    fn from_iter<I: IntoIterator<Item = Monomer>>(iter: I) -> Self {
        let mut k = Self::new();

        for i in iter {
            k.add(i.into());
        }
        k
    }
}

impl Monomer {
    fn parse(m: u8) -> Result<Self, ValidityError> {
        match m {
            b'A' => Ok(Self::A(m)),
            b'C' => Ok(Self::C(m)),
            b'G' => Ok(Self::G(m)),
            b'T' => Ok(Self::T(m)),
            _ => Err(ValidityError::InvalidByte),
        }
    }

    fn complement(self) -> Self {
        match self {
            Self::A(_) => Self::from(b'T'),
            Self::C(_) => Self::from(b'G'),
            Self::G(_) => Self::from(b'C'),
            Self::T(_) => Self::from(b'A'),
        }
    }

    fn pack(self) -> u64 {
        match self {
            Self::A(_) => 0,
            Self::C(_) => 1,
            Self::G(_) => 2,
            Self::T(_) => 3,
        }
    }

    fn unpack(bit: u64) -> Self {
        match bit {
            0 => Self::from(b'A'),
            1 => Self::from(b'C'),
            2 => Self::from(b'G'),
            _ => Self::from(b'T'),
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

    fn pack(&mut self, elem: u8) {
        self.shift();
        let mask = Monomer::pack(elem.into());
        self.0 |= mask
    }

    fn shift(&mut self) {
        self.0 <<= 2
    }
}

impl FromIterator<u8> for Bitpack {
    fn from_iter<I: IntoIterator<Item = u8>>(iter: I) -> Self {
        let mut c = Self::new();

        for i in iter {
            c.pack(i)
        }
        c
    }
}

/// Unpack bitpacked k-mer data
#[derive(Hash, PartialEq, Eq)]
pub struct Unpack(pub Vec<u8>);

impl Unpack {
    fn new() -> Self {
        Unpack(Vec::new())
    }

    fn add(&mut self, elem: u8) {
        self.0.push(elem);
    }

    pub(crate) fn bit(bit: u64, k: usize) -> Self {
        (0..k)
            .into_iter()
            .map(|i| bit.isolate(i, k))
            .map(|bit| bit.replace())
            .map(Monomer::unpack)
            .map(|monomer| monomer.into())
            .collect()
    }
}

impl FromIterator<u8> for Unpack {
    fn from_iter<I: IntoIterator<Item = u8>>(iter: I) -> Self {
        let mut c = Unpack::new();

        for i in iter {
            c.add(i)
        }
        c
    }
}

/// Convert a DNA string slice into its [reverse compliment](https://en.wikipedia.org/wiki/Complementarity_(molecular_biology)#DNA_and_RNA_base_pair_complementarity).
pub(crate) struct RevComp(pub Vec<u8>);

impl FromIterator<u8> for RevComp {
    fn from_iter<I: IntoIterator<Item = u8>>(iter: I) -> Self {
        let mut c = RevComp::new();

        for i in iter {
            c.add(i)
        }
        c
    }
}

impl RevComp {
    fn new() -> Self {
        Self(vec![])
    }

    fn add(&mut self, elem: u8) {
        self.0.push(elem)
    }

    pub(crate) fn from_kmer(kmer: &Kmer) -> Self {
        kmer.0
            .iter()
            .rev()
            .map(|byte| Monomer::from(*byte).complement().into())
            .collect()
    }
}

#[cfg(test)]
pub mod test {
    use super::*;

    #[test]
    fn test_from_substring() {
        let sub = &[b'G', b'A', b'T', b'T', b'A', b'C', b'A'];
        let k = Kmer::from_sub(sub).unwrap();
        insta::assert_snapshot!(format!("{:?}", k), @"Kmer([71, 65, 84, 84, 65, 67, 65])");
    }

    #[test]
    fn test_parse_valid_byte() {
        let b = b'N';
        assert!(Monomer::parse(b).is_err());
    }

    #[test]
    fn from_substring_returns_err_for_invalid_substring() {
        let sub = &[b'N'];
        let k = Kmer::from_sub(sub);
        assert!(k.is_err());
    }

    #[test]
    fn find_invalid_works() {
        let dna = "NACNN".as_bytes();
        let ans = Kmer::find_invalid(dna);
        assert_eq!(4, ans);
        assert_eq!(&b'N', dna.iter().collect::<Vec<_>>()[ans]);

        let dna = "NACNG".as_bytes();
        let ans = Kmer::find_invalid(dna);
        assert_eq!(3, ans);
        assert_eq!(&b'N', dna.iter().collect::<Vec<_>>()[ans]);

        let dna = "NANTG".as_bytes();
        let ans = Kmer::find_invalid(dna);
        assert_eq!(2, ans);
        assert_eq!(&b'N', dna.iter().collect::<Vec<_>>()[ans]);

        let dna = "NNCTG".as_bytes();
        let ans = Kmer::find_invalid(dna);
        assert_eq!(1, ans);
        assert_eq!(&b'N', dna.iter().collect::<Vec<_>>()[ans]);

        let dna = "NACTG".as_bytes();
        let ans = Kmer::find_invalid(dna);
        assert_eq!(0, ans);
        assert_eq!(&b'N', dna.iter().collect::<Vec<_>>()[ans]);
    }
}
