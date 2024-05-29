use ark_ff::{BigInteger, PrimeField};
use rand::Rng;
use sha3::{Digest, Keccak256};

use polynomials::{
    multilinear_polynomial::traits::MultilinearPolynomialTrait,
    univariate_polynomial::UnivariatePolynomial,
};

pub struct Transcript {
    hasher: Keccak256,
}

impl Transcript {
    pub fn new() -> Self {
        Self {
            hasher: Keccak256::new(),
        }
    }

    pub fn append(&mut self, new_data: &[u8]) {
        self.hasher.update(&mut new_data.to_owned());
    }

    pub fn sample_challenge(&mut self) -> [u8; 32] {
        let mut result_hash = [0_u8; 32];
        result_hash.copy_from_slice(&self.hasher.finalize_reset());
        result_hash.reverse();
        self.hasher.update(result_hash);
        result_hash
    }

    pub fn sample_field_element<F: PrimeField>(&mut self) -> F {
        let challenge = self.sample_challenge();
        F::from_be_bytes_mod_order(&challenge)
    }

    pub fn add_univariate_poly<F: PrimeField>(
        &mut self,
        univariate_poly: &UnivariatePolynomial<F>,
    ) {
        for i in 0..univariate_poly.coefficients.len() {
            self.append(
                univariate_poly.coefficients[i]
                    .into_bigint()
                    .to_bytes_be()
                    .as_slice(),
            )
        }
    }

    pub fn add_multivariate_poly<F: PrimeField, MPT: MultilinearPolynomialTrait<F>>(
        &mut self,
        multilinear_poly: &MPT,
    ) {
        self.append(&multilinear_poly.to_bytes());
    }

    pub fn sample_n_field_elements<F: PrimeField>(&mut self, n: usize) -> Vec<F> {
        let mut res = vec![];

        for _ in 0..n {
            let r_i = self.sample_field_element();
            res.push(r_i);
        }

        res
    }

    pub fn sample_element_in_range_n(&self, n: u64) -> u64 {
        let mut rng = rand::thread_rng();
        rng.gen_range(0..n)
    }
}
