/// Creating a valid k-mer bytestring.
#[derive(Debug, PartialEq)]
pub struct Kmer(pub Vec<u8>);

impl Kmer {
    pub fn new(sub: &[u8]) -> Option<Kmer> {
        if !sub.contains(&b'N') {
            let valid_kmer = sub.to_vec();
            Some(Kmer(valid_kmer))
        } else {
            None
        }
    }

    /// Find the index of the rightmost invalid byte in an invalid bytestring.
    pub fn find_invalid(sub: &[u8]) -> usize {
        match sub
            .iter()
            .rposition(|byte| ![b'A', b'C', b'G', b'T'].contains(byte))
        {
            Some(rightmost_invalid_byte_index) => rightmost_invalid_byte_index,
            None => panic!("Valid bytestring passed to `find_invalid`, which is a bug."),
        }
    }
}
