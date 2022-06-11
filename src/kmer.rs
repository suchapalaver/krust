/// Creating a valid k-mer bytestring.
#[derive(Debug, PartialEq)]
pub struct Kmer(pub Vec<u8>);

impl Kmer {
    pub fn new() -> Kmer {
        Kmer(Vec::new())
    }

    fn add(&mut self, elem: u8) {
        self.0.push(elem)
    }

    pub fn from_substring(sub: &[u8]) -> Result<Self, ()> {
        sub.iter()
            .map(|b| b.parse_valid_byte())
            .collect()
    }

    /// Find the index of the rightmost invalid byte in an invalid bytestring.
    pub fn find_invalid(sub: &[u8]) -> usize {
        match sub
            .iter()
	    .rposition(|byte| byte.parse_valid_byte().is_err())
        {
            Some(rightmost_invalid_byte_index) => rightmost_invalid_byte_index,
            None => panic!("Valid bytestring passed to `find_invalid`, which is a bug."),
        }
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
    fn parse_valid_byte(self) -> Result<Self, ()>
    where
        Self: Sized;
}

impl Validity for u8 {
    fn parse_valid_byte(self) -> Result<Self, ()> {
        if [b'A', b'C', b'G', b'T'].contains(&self) {
            Ok(self)
	} else {
            Err(())
        }
    }
}

#[cfg(test)]
pub mod test {
    use super::*;

    #[test]
    fn test_from_substring() {
        let sub = &[b'C', b'A', b'G', b'T'];
        match Kmer::from_substring(sub) {
            Ok(k) => insta::assert_snapshot!(format!("{:?}", k), @"Kmer([67, 65, 71, 84])"),
            Err(()) => panic!("this should not happen"),
        }
    }

     #[test]
    fn test_parse_valid_byte() {
	let sub = &[b'C', b'N', b'G', b'T'];
	assert!(sub[1].parse_valid_byte().is_err());
    }

    #[test]
    fn test_from_substring_returns_err_for_invalid_substring() {
	let sub = &[b'C', b'N', b'G', b'T'];
	let k = Kmer::from_substring(sub);
	assert!(k.is_err());
    }
}
