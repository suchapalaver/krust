use krust::startup::Kmer;

#[test]
fn find_invalid_works1() {
    let dna = "NACNN".as_bytes();
    let ans = Kmer::find_invalid(dna);
    assert_eq!(4, ans);
    assert_eq!(&b'N', dna.iter().collect::<Vec<_>>()[ans]);
}

#[test]
fn find_invalid_works2() {
    let dna = "NACNG".as_bytes();
    let ans = Kmer::find_invalid(dna);
    assert_eq!(3, ans);
    assert_eq!(&b'N', dna.iter().collect::<Vec<_>>()[ans]);
}

#[test]
fn find_invalid_works3() {
    let dna = "NANTG".as_bytes();
    let ans = Kmer::find_invalid(dna);
    assert_eq!(2, ans);
    assert_eq!(&b'N', dna.iter().collect::<Vec<_>>()[ans]);
}

#[test]
fn find_invalid_works4() {
    let dna = "NNCTG".as_bytes();
    let ans = Kmer::find_invalid(dna);
    assert_eq!(1, ans);
    assert_eq!(&b'N', dna.iter().collect::<Vec<_>>()[ans]);
}

#[test]
fn find_invalid_works5() {
    let dna = "NACTG".as_bytes();
    let ans = Kmer::find_invalid(dna);
    assert_eq!(0, ans);
    assert_eq!(&b'N', dna.iter().collect::<Vec<_>>()[ans]);
}

#[test]
#[should_panic(expected = "Valid bytestring passed to `find_invalid`, which is a bug.")]
fn find_invalid_panics_when_passed_valid_kmer() {
    let dna = "CACTG".as_bytes();
    let _ans = Kmer::find_invalid(dna);
}
