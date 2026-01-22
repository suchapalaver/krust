use bytes::Bytes;
use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};
use kmerust::kmer::{unpack_to_bytes, Kmer, KmerLength};
use kmerust::run::count_kmers;
use std::io::Write;
use tempfile::NamedTempFile;

fn bench_from_sub(c: &mut Criterion) {
    let mut group = c.benchmark_group("Kmer::from_sub");

    for k in [5, 11, 21, 31] {
        let seq = "ACGT".repeat(k / 4 + 1);
        let bytes = Bytes::copy_from_slice(&seq.as_bytes()[..k]);

        group.bench_with_input(BenchmarkId::from_parameter(k), &bytes, |b, bytes| {
            b.iter(|| Kmer::from_sub(black_box(bytes.clone())))
        });
    }

    group.finish();
}

fn bench_pack(c: &mut Criterion) {
    let mut group = c.benchmark_group("Kmer::pack");

    for k in [5, 11, 21, 31] {
        let seq = "ACGT".repeat(k / 4 + 1);
        let bytes = Bytes::copy_from_slice(&seq.as_bytes()[..k]);

        group.bench_with_input(BenchmarkId::from_parameter(k), &bytes, |b, bytes| {
            b.iter(|| {
                let kmer = Kmer::from_sub(bytes.clone()).unwrap();
                black_box(kmer.pack())
            })
        });
    }

    group.finish();
}

fn bench_canonical(c: &mut Criterion) {
    let mut group = c.benchmark_group("Kmer::canonical");

    for k in [5, 11, 21, 31] {
        let seq = "ACGT".repeat(k / 4 + 1);
        let bytes = Bytes::copy_from_slice(&seq.as_bytes()[..k]);

        group.bench_with_input(BenchmarkId::from_parameter(k), &bytes, |b, bytes| {
            b.iter(|| {
                let kmer = Kmer::from_sub(bytes.clone()).unwrap();
                black_box(kmer.pack().canonical())
            })
        });
    }

    group.finish();
}

fn bench_unpack(c: &mut Criterion) {
    let mut group = c.benchmark_group("unpack_to_bytes");

    for k in [5, 11, 21, 31] {
        let seq = "ACGT".repeat(k / 4 + 1);
        let bytes = Bytes::copy_from_slice(&seq.as_bytes()[..k]);
        let packed = Kmer::from_sub(bytes).unwrap().pack();
        let packed_bits = packed.packed_bits();
        let k_len = KmerLength::new(k).unwrap();

        group.bench_with_input(
            BenchmarkId::from_parameter(k),
            &(packed_bits, k_len),
            |b, &(bits, k_len)| b.iter(|| black_box(unpack_to_bytes(bits, k_len))),
        );
    }

    group.finish();
}

fn bench_count_kmers_small(c: &mut Criterion) {
    let mut group = c.benchmark_group("count_kmers");

    // Create a small test file
    let mut file = NamedTempFile::new().unwrap();
    for i in 0..100 {
        writeln!(file, ">seq{i}").unwrap();
        writeln!(file, "{}", "ACGTACGTACGTACGTACGTACGTACGTACGT".repeat(10)).unwrap();
    }
    let path = file.path().to_path_buf();

    for k in [5, 11, 21] {
        group.bench_with_input(BenchmarkId::from_parameter(k), &k, |b, &k| {
            b.iter(|| count_kmers(black_box(&path), black_box(k)))
        });
    }

    group.finish();
}

criterion_group!(
    benches,
    bench_from_sub,
    bench_pack,
    bench_canonical,
    bench_unpack,
    bench_count_kmers_small,
);

criterion_main!(benches);
