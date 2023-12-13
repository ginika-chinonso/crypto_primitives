use ark_ff::{BigInteger, PrimeField};
use sha3::{Digest, Keccak256};

use crate::polynomials::{
    multilinear_poly::MultilinearPolynomial, univariate_poly::UnivariatePolynomial,
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

    pub fn add_multivariate_poly<F: PrimeField>(
        &mut self,
        multilinear_poly: &MultilinearPolynomial<F>,
    ) {
        for i in 0..multilinear_poly.terms.len() {
            self.append(
                multilinear_poly.terms[i]
                    .coefficient
                    .into_bigint()
                    .to_bytes_be()
                    .as_slice(),
            );
            self.append(
                multilinear_poly.terms[i]
                    .vars
                    .iter()
                    .map(|a| *a as u8)
                    .collect::<Vec<u8>>()
                    .as_slice(),
            );
        }
    }
}
