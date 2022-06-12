custom_error::custom_error! { pub ValidityError
    InvalidByte = "",
}

/// Creating a valid k-mer bytestring.
#[derive(Debug, PartialEq)]
pub struct Kmer(pub Vec<u8>);

/// Find the canonical kmer
/// --the alphabetically smaller of the substring and its reverse complement.
pub struct CanonicalKmer(pub Vec<u8>);

impl Kmer {
    pub fn new() -> Kmer {
        Kmer(Vec::new())
    }

    fn add(&mut self, elem: u8) {
        self.0.push(elem)
    }

    pub fn from_substring(sub: &[u8]) -> Result<Kmer, ValidityError> {
        sub.iter().map(|b| b.parse_valid_byte()).collect()
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
            .rposition(|byte| byte.parse_valid_byte().is_err())
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
            k.add(i);
        }
        k
    }
}

trait Validity {
    fn parse_valid_byte(self) -> Result<Self, ValidityError>
    where
        Self: Sized;
}

impl Validity for u8 {
    fn parse_valid_byte(self) -> Result<Self, ValidityError> {
        if [b'A', b'C', b'G', b'T'].contains(&self) {
            Ok(self)
        } else {
            Err(ValidityError::InvalidByte)
        }
    }
}

#[cfg(test)]
pub mod test {
    use super::*;

    #[test]
    fn test_from_substring() {
        let sub = &[b'C', b'A', b'G', b'T', b'G'];
        let k = Kmer::from_substring(sub).unwrap();
        insta::assert_snapshot!(format!("{:?}", k), @"Kmer([67, 65, 71, 84, 71])");
    }

    #[test]
    fn test_parse_valid_byte() {
        let sub = &[b'C', b'N', b'G', b'T', b'G'];
        assert!(sub[1].parse_valid_byte().is_err());
    }

    #[test]
    fn test_from_substring_returns_err_for_invalid_substring() {
        let sub = &[b'C', b'N', b'G', b'T'];
        let k = Kmer::from_substring(sub);
        assert!(k.is_err());
    }
}
