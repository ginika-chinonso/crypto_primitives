use std::{error::Error, marker::PhantomData};

use ark_ff::PrimeField;
use fiat_shamir_transcript::Transcript;
use polynomials::multilinear_polynomial::traits::MultilinearPolynomialTrait;

use super::prover::{SumcheckProductProof, SumcheckProof};

pub struct Verifier<F> {
    pub _marker: PhantomData<F>,
}

impl<F: PrimeField + From<i32>> Verifier<F> {
    pub fn verify_sumcheck(proof: SumcheckProof<F>) -> Result<bool, Box<dyn Error>> {
        let mut transcript = Transcript::new();
        transcript.append(&proof.init_poly.get(0).unwrap().to_bytes());

        let mut challenges: Vec<(usize, F)> = vec![];

        let mut sum = proof.claimed_sum;

        for i in 0..proof.round_polys.len() {
            if sum
                != proof.round_polys[i].evaluate(&vec![(1, F::zero())])
                    + proof.round_polys[i].evaluate(&vec![(1, F::one())])
            {
                return Err(format!("Verifier check for round {} failed", i).into());
            }

            let challenge = transcript.sample_field_element();

            challenges.push((i + 1, challenge));

            sum = proof.round_polys[i].evaluate(&vec![(1, challenge)]);

            transcript.append(&proof.round_polys[i].to_bytes());
        }

        if proof.init_poly.get(0).unwrap().evaluate(&challenges) != sum {
            return Err(format!("Verifier oracle query failed").into());
        }

        return Ok(true);
    }

    pub fn verify_sumcheck_product(proof: SumcheckProductProof<F>) -> Result<bool, Box<dyn Error>> {
        let (f_poly, g_poly) = (&proof.init_poly[0], &proof.init_poly[1]);

        let mut transcript = Transcript::new();
        transcript.append(&f_poly.to_bytes());
        transcript.append(&g_poly.to_bytes());

        let mut claimed_sum = proof.claimed_sum;

        let mut challenges = vec![];

        for i in 0..f_poly.num_of_vars {
            let round_poly = &proof.round_polys[i];

            dbg!(&claimed_sum);
            dbg!(&proof.round_polys[i]);
            dbg!(round_poly.evaluate(F::zero()) + round_poly.evaluate(F::one()));

            if claimed_sum != (round_poly.evaluate(F::zero()) + round_poly.evaluate(F::one())) {
                return Err(format!("Verifier check for round {} failed", i).into());
            }

            let challenge = transcript.sample_field_element();

            dbg!("i: {}, challenge: {}", i, challenge);

            claimed_sum = round_poly.evaluate(challenge);

            transcript.append(&round_poly.to_bytes());
            challenges.push((i + 1, challenge));
        }

        if claimed_sum != (f_poly.evaluate(&challenges) * g_poly.evaluate(&challenges)) {
            return Err(format!("Verifier oracle check failed").into());
        }

        Ok(true)
    }
}
