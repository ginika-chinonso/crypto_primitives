use crate::polynomials::univariate_poly::UnivariatePolynomial;
use ark_ff::PrimeField;

// Checks if two vectors are the same using reed solomon fingerprinting
pub fn reed_solomon_fingerprinting<F: PrimeField + std::convert::From<i32>>(
    a: Vec<F>,
    b: Vec<F>,
    challenge: F,
) -> bool {
    let a_poly = UnivariatePolynomial::new(a);
    let b_poly = UnivariatePolynomial::new(b);

    a_poly.evaluate(challenge) == b_poly.evaluate(challenge)
}

#[cfg(test)]
mod tests {
    use super::reed_solomon_fingerprinting;
    use ark_ff::{Fp64, MontBackend, MontConfig};
    use ark_std::UniformRand;

    #[derive(MontConfig)]
    #[modulus = "17"]
    #[generator = "3"]
    pub struct FqConfig;
    pub type Fq = Fp64<MontBackend<FqConfig, 1>>;

    #[test]
    fn test_reed_solomon_same() {
        let mut rng = ark_std::test_rng();
        let poly_a = vec![Fq::rand(&mut rng), Fq::rand(&mut rng), Fq::rand(&mut rng)];
        let challenge = Fq::rand(&mut rng);

        let res = reed_solomon_fingerprinting(poly_a.clone(), poly_a, challenge);
        assert!(res);
    }

    #[test]
    fn test_reed_solomon_diff() {
        let mut rng = ark_std::test_rng();
        let poly_a = vec![Fq::rand(&mut rng), Fq::rand(&mut rng), Fq::rand(&mut rng)];
        let poly_b = vec![Fq::rand(&mut rng), Fq::rand(&mut rng), Fq::rand(&mut rng)];
        let challenge = Fq::rand(&mut rng);

        let res = reed_solomon_fingerprinting(poly_a, poly_b, challenge);
        assert!(!res);
    }
}
