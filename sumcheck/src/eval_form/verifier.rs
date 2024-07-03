use std::{error::Error, marker::PhantomData};

use ark_ff::PrimeField;
use fiat_shamir_transcript::Transcript;
use polynomials::multilinear_polynomial::traits::MultilinearPolynomialTrait;

use super::prover::SumcheckProof;

pub struct Verifier<F> {
    pub _marker: PhantomData<F>
}

impl <F: PrimeField>Verifier<F> {

    pub fn verify(proof: SumcheckProof<F>) -> Result<bool, Box<dyn Error>> {
        
        let mut transcript = Transcript::new();
        transcript.append(&proof.init_poly.to_bytes());

        let mut challenges: Vec<(usize, F)> = vec![];

        let mut sum = proof.claimed_sum;

        for i in 0..proof.round_polys.len() {

            if sum != proof.round_polys[i].evaluate(&vec![(1, F::zero())]) + proof.round_polys[i].evaluate(&vec![(1, F::one())]) {
                return Err(format!("Verifier check for round {} failed", i).into());
            }

            let challenge = transcript.sample_field_element();
            
            challenges.push((i + 1, challenge));
            
            sum = proof.round_polys[i].evaluate(&vec![(1, challenge)]);
            
            transcript.append(&proof.round_polys[i].to_bytes());
        }

        if proof.init_poly.evaluate(&challenges) != sum {
            return Err(format!("Verifier oracle query failed").into());
        }

        return Ok(true);
    }
}