use bitgauss::{bitmatrix::BitMatrix, bitvec::*};
use criterion::{criterion_group, criterion_main, BatchSize, Criterion};
use rand::{rngs::SmallRng, Rng, SeedableRng};
// use std::hint::black_box;

fn bitvec_ops(c: &mut Criterion) {
    let mut group = c.benchmark_group("bitvec_ops");
    let sz = 1563; // ~100K bits
    let mut rng = SmallRng::seed_from_u64(1);
    let vec1: BitVec = (0..sz).map(|_| rng.random::<BitBlock>()).collect();
    let vec2: BitVec = (0..sz).map(|_| rng.random::<BitBlock>()).collect();
    let vec3: BitVec = (0..(2 * sz)).map(|_| rng.random::<BitBlock>()).collect();

    group.bench_function("bitvec_xor", |b| {
        b.iter_batched_ref(
            || vec1.clone(),
            |bv| {
                bv.bitxor_assign(&vec2);
            },
            BatchSize::LargeInput,
        )
    });

    group.bench_function("bitvec_xor_slice", |b| {
        b.iter_batched_ref(
            || vec3.clone(),
            |bv| bv.xor_range(10, sz + 10, sz - 20),
            BatchSize::LargeInput,
        )
    });

    group.bench_function("bitvec_and", |b| {
        b.iter_batched_ref(
            || vec1.clone(),
            |bv| {
                bv.bitand_assign(&vec2);
            },
            BatchSize::LargeInput,
        )
    });

    group.bench_function("bitvec_dot", |b| b.iter(|| vec1.dot(&vec2)));
}

fn gauss(c: &mut Criterion) {
    let mut group = c.benchmark_group("gauss");
    let mut rng = SmallRng::seed_from_u64(1);
    let m = BitMatrix::random(&mut rng, 1000, 1000);
    group.bench_function("gauss_utri", |b| {
        b.iter_batched_ref(|| m.clone(), |m| m.gauss(false), BatchSize::LargeInput)
    });

    group.bench_function("gauss_full", |b| {
        b.iter_batched_ref(|| m.clone(), |m| m.gauss(true), BatchSize::LargeInput)
    });
}

fn big_gauss(c: &mut Criterion) {
    let mut group = c.benchmark_group("big_gauss");
    group.sample_size(10);
    let mut rng = SmallRng::seed_from_u64(1);
    let m = BitMatrix::random(&mut rng, 10_000, 10_000);
    group.bench_function("big_gauss_utri", |b| {
        b.iter_batched_ref(|| m.clone(), |m| m.gauss(false), BatchSize::LargeInput)
    });

    group.bench_function("big_gauss_full", |b| {
        b.iter_batched_ref(|| m.clone(), |m| m.gauss(true), BatchSize::LargeInput)
    });
}

fn transpose(c: &mut Criterion) {
    let mut group = c.benchmark_group("transpose");
    let mut rng = SmallRng::seed_from_u64(1);
    let m = BitMatrix::random(&mut rng, 1000, 1000);

    group.bench_function("naive", |b| {
        b.iter_batched_ref(
            || m.clone(),
            |m| {
                for i in 0..m.rows() {
                    for j in 0..m.cols() {
                        let b0 = m.bit(i, j);
                        m.set_bit(i, j, m.bit(j, i));
                        m.set_bit(j, i, b0);
                    }
                }
            },
            BatchSize::LargeInput,
        )
    });

    group.bench_function("inplace", |b| {
        b.iter_batched_ref(
            || m.clone(),
            |m| m.transpose_inplace(),
            BatchSize::LargeInput,
        )
    });

    group.bench_function("copied", |b| {
        b.iter_batched_ref(|| m.clone(), |m| m.transposed(), BatchSize::LargeInput)
    });
}

fn mult(c: &mut Criterion) {
    let mut group = c.benchmark_group("matrix_mult");
    let mut rng = SmallRng::seed_from_u64(1);
    let m1 = BitMatrix::random(&mut rng, 200, 200);
    let m2 = BitMatrix::random(&mut rng, 200, 200);
    group.bench_function("medium", |b| b.iter(|| &m1 * &m2));

    let m1 = BitMatrix::random(&mut rng, 1000, 1000);
    let m2 = BitMatrix::random(&mut rng, 1000, 1000);
    group.bench_function("large", |b| b.iter(|| &m1 * &m2));

    group.sample_size(10);

    let m1 = BitMatrix::random(&mut rng, 10_000, 10_000);
    let m2 = BitMatrix::random(&mut rng, 10_000, 10_000);
    group.bench_function("huge", |b| b.iter(|| &m1 * &m2));
}

criterion_group!(benches, bitvec_ops, gauss, big_gauss, transpose, mult);
criterion_main!(benches);
