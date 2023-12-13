use ark_poly::{DenseMultilinearExtension, MultilinearExtension};
use criterion::{criterion_group, criterion_main, Criterion};
use crypto_primitives::polynomials::multilinear_poly::{
    MultilinearMonomial, MultilinearPolynomial,
};

use std::vec;

use ark_ff::{Fp64, MontBackend, MontConfig};

#[derive(MontConfig)]
#[modulus = "17"]
#[generator = "3"]
pub struct FqConfig;
pub type Fq = Fp64<MontBackend<FqConfig, 1>>;

fn multilinear_poly_benchmark(c: &mut Criterion) {
    let poly = MultilinearPolynomial::new(vec![
        MultilinearMonomial::new(Fq::from(3), vec![true, true, false]),
        MultilinearMonomial::new(Fq::from(4), vec![false, true, true]),
    ]);
    let poly1 = MultilinearPolynomial::new(vec![
        MultilinearMonomial::new(Fq::from(5), vec![false, true, false]),
        MultilinearMonomial::new(Fq::from(9), vec![true, false, true]),
    ]);

    let poly3 = MultilinearPolynomial::<Fq>::interpolate(vec![
        Fq::from(2),
        Fq::from(4),
        Fq::from(6),
        Fq::from(8),
        Fq::from(10),
        Fq::from(15),
        Fq::from(3),
        Fq::from(1),
    ]);

    c.bench_function("multilinear evaluation", |b| {
        b.iter(|| {
            poly3
                .clone()
                .evaluate(vec![(0, Fq::from(3)), (1, Fq::from(9)), (2, Fq::from(8))])
        })
    });

    c.bench_function("multilinear addition", |b| {
        b.iter(|| poly.clone().add(poly1.clone()))
    });

    c.bench_function("multilinear interpolate", |b| {
        b.iter(|| {
            MultilinearPolynomial::<Fq>::interpolate(vec![
                Fq::from(2),
                Fq::from(4),
                Fq::from(6),
                Fq::from(8),
                Fq::from(10),
            ])
        })
    });

    c.bench_function("multilinear multiplication", |b| {
        b.iter(|| poly.clone().multiply(&mut poly1.clone()))
    });

    c.bench_function("multilinear scalar mul", |b| {
        b.iter(|| poly.clone().scalar_mul(Fq::from(8)))
    });

    c.bench_function("arkworks evaluate", |b| {
        b.iter(|| {
            DenseMultilinearExtension::from_evaluations_vec(
                3,
                vec![
                    Fq::from(2),
                    Fq::from(4),
                    Fq::from(6),
                    Fq::from(8),
                    Fq::from(10),
                    Fq::from(15),
                    Fq::from(3),
                    Fq::from(1),
                ],
            )
            .evaluate(&vec![Fq::from(3), Fq::from(9), Fq::from(8)])
        })
    });
}

criterion_group!(multilinear_poly, multilinear_poly_benchmark);
criterion_main!(multilinear_poly);
