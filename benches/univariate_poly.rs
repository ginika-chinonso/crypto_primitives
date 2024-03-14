use criterion::{criterion_group, criterion_main, Criterion};
use polynomials::univariate_polynomial::UnivariatePolynomial;

use std::vec;

use ark_ff::{Fp64, MontBackend, MontConfig};

#[derive(MontConfig)]
#[modulus = "17"]
#[generator = "3"]
pub struct FqConfig;
pub type Fq = Fp64<MontBackend<FqConfig, 1>>;

fn univariate_poly_benchmark(c: &mut Criterion) {
    let new_polynomial = UnivariatePolynomial::new(vec![
        Fq::from(134546355),
        Fq::from(2444543),
        Fq::from(3464556),
    ]);
    let poly1 = UnivariatePolynomial::new(vec![Fq::from(1000), Fq::from(2000), Fq::from(3005)]);
    let poly2 = UnivariatePolynomial::new(vec![
        Fq::from(1000000),
        Fq::from(2000000),
        Fq::from(308677675),
    ]);

    c.bench_function("evaluate poly", |b| {
        b.iter(|| new_polynomial.evaluate(Fq::from(6)))
    });

    c.bench_function("add poly", |b| b.iter(|| poly1.clone() + poly2.clone()));

    c.bench_function("mul poly", |b| b.iter(|| poly1.clone() * poly2.clone()));

    c.bench_function("interpolate poly", |b| {
        b.iter(|| {
            UnivariatePolynomial::interpolate(
                poly1.coefficients.clone(),
                poly2.coefficients.clone(),
            )
        })
    });
}

criterion_group!(univariate_poly, univariate_poly_benchmark);
criterion_main!(univariate_poly);
