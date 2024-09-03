use ark_ff::{FftField, PrimeField};
use polynomials::univariate_polynomial::UnivariatePolynomial;

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

// Reed Solomon code word
// Generates a reed solomon code word for a message
pub fn generate_codeword<F: PrimeField + FftField>(
    message: &Vec<F>,
    blow_up_factor: usize,
) -> Vec<F> {
    let domain_generator: F = FftField::get_root_of_unity(message.len() as u64).unwrap();
    let domain = (0..message.len())
        .into_iter()
        .map(|val| domain_generator.pow([val as u64]))
        .collect();
    let poly = UnivariatePolynomial::interpolate(&domain, &message.to_vec());
    let blow_up_domain_generator: F =
        FftField::get_root_of_unity((message.len() * blow_up_factor) as u64).unwrap();
    let blow_up_domain: Vec<F> = (0..(message.len() * blow_up_factor))
        .into_iter()
        .map(|val| blow_up_domain_generator.pow([val as u64]))
        .collect();
    let res = blow_up_domain
        .into_iter()
        .map(|val| poly.evaluate(val))
        .collect();
    res
}

#[cfg(test)]
mod tests {
    use super::reed_solomon_fingerprinting;
    use ark_bn254::Fq;
    use ark_std::UniformRand;

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

    #[test]
    pub fn test_generate_codeword() {
        todo!()
    }
}
