/// Find the canonical kmer
/// --the alphabetically smaller of the substring and its reverse complement.
pub struct CanonicalKmer(pub Vec<u8>);

impl From<(Vec<u8>, Vec<u8>)> for CanonicalKmer {
    fn from(comp: (Vec<u8>, Vec<u8>)) -> Self {
        let (reverse_complement, kmer) = (comp.0, comp.1);
        let canonical_kmer = if reverse_complement.cmp(&kmer) == std::cmp::Ordering::Less {
            reverse_complement
        } else {
            kmer
        };
        CanonicalKmer(canonical_kmer)
    }
}
