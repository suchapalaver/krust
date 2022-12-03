custom_error::custom_error! { pub ValidityError
    InvalidByte = "Not a valid monomer",
}

/// A valid k-mer
#[derive(Debug, Eq, PartialEq)]
pub struct Kmer(pub Vec<u8>);

/// Find the canonical kmer
///
/// # Notes 
/// The alphabetically smaller of the substring and its reverse complement
pub struct CanonicalKmer(pub Vec<u8>);

impl Kmer {
    pub fn new() -> Kmer {
        Kmer(Vec::new())
    }

    fn add(&mut self, elem: Monomer) {
        self.0.push(elem.into())
    }

    pub fn from_substring(sub: &[u8]) -> Result<Kmer, ValidityError> {
        sub.iter().map(|b| parse_valid_byte(*b)).collect()
    }

    pub fn get_canonical_kmer(reverse_complement: Vec<u8>, kmer: Vec<u8>) -> Self {
        let rc = {
            if reverse_complement.cmp(&kmer) == std::cmp::Ordering::Less {
                reverse_complement
            } else {
                kmer
            }
        };
        rc.into_iter().collect()
    }

    pub fn find_invalid_byte_index(sub: &[u8]) -> usize {
        sub.iter()
            .rposition(|byte| parse_valid_byte(*byte).is_err())
            .unwrap()
    }
}

impl Default for Kmer {
    fn default() -> Self {
        Self::new()
    }
}

impl FromIterator<u8> for Kmer {
    fn from_iter<I: IntoIterator<Item = u8>>(iter: I) -> Self {
        let mut k = Kmer::new();

        for i in iter {
            k.add(i.into());
        }
        k
    }
}

enum Monomer {
    A(u8), 
    C(u8), 
    G(u8), 
    T(u8)
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
        let mut k = Kmer::new();

        for i in iter {
            k.add(i);
        }
        k
    }
}

fn parse_valid_byte(m: u8) -> Result<Monomer, ValidityError> {
    match m {
        b'A' => Ok(Monomer::A(m)),
        b'C' => Ok(Monomer::C(m)),
        b'G' => Ok(Monomer::G(m)),
        b'T' => Ok(Monomer::T(m)),
        _ => Err(ValidityError::InvalidByte),
    }
}

#[cfg(test)]
pub mod test {
    use super::*;

    #[test]
    fn test_from_substring() {
        let sub = &[b'C', b'A', b'G', b'T'];
        let k = Kmer::from_substring(sub).unwrap();
        insta::assert_snapshot!(format!("{:?}", k), @"Kmer([67, 65, 71, 84])");
    }

    #[test]
    fn test_parse_valid_byte() {
        let b = b'N';
        assert!(parse_valid_byte(b).is_err());
    }

    #[test]
    fn test_from_substring_returns_err_for_invalid_substring() {
        let sub = &[b'N'];
        let k = Kmer::from_substring(sub);
        assert!(k.is_err());
    }
}
