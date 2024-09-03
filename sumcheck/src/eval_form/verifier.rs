use std::{error::Error, marker::PhantomData};

use ark_ff::PrimeField;
use fiat_shamir_transcript::Transcript;
use polynomials::multilinear_polynomial::{eval_form::MLE, traits::MultilinearPolynomialTrait};

use crate::universal_mle::universal_mle::SumcheckProof;

pub struct Verifier<F> {
    pub _marker: PhantomData<F>,
}

impl<F: PrimeField + From<i32>> Verifier<F> {
    pub fn verify_sumcheck(
        init_poly: MLE<F>,
        proof: SumcheckProof<F>,
    ) -> Result<bool, Box<dyn Error>> {
        let mut transcript = Transcript::new();
        transcript.append(&init_poly.to_bytes());

        let mut challenges: Vec<(usize, F)> = vec![];

        let mut sum = proof.claimed_sum;

        for i in 0..proof.round_polys.len() {
            if sum
                != proof.round_polys[i].evaluate(F::zero())
                    + proof.round_polys[i].evaluate(F::one())
            {
                return Err(format!("Verifier check for round {} failed", i).into());
            }

            let challenge = transcript.sample_field_element();

            challenges.push((i + 1, challenge));

            sum = proof.round_polys[i].evaluate(challenge);

            transcript.append(&proof.round_polys[i].to_bytes());
        }

        if init_poly.evaluate(&challenges) != sum {
            return Err(format!("Verifier oracle query failed").into());
        }

        return Ok(true);
    }

    pub fn verify_sumcheck_product(
        init_poly: Vec<MLE<F>>,
        proof: SumcheckProof<F>,
    ) -> Result<bool, Box<dyn Error>> {
        let mut transcript = Transcript::new();
        transcript.append(&init_poly[0].to_bytes());
        transcript.append(&init_poly[1].to_bytes());

        let mut claimed_sum = proof.claimed_sum;

        let mut challenges = vec![];

        for i in 0..init_poly[0].num_of_vars {
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

        if claimed_sum != (init_poly[0].evaluate(&challenges) * init_poly[1].evaluate(&challenges))
        {
            return Err(format!("Verifier oracle check failed").into());
        }

        Ok(true)
    }
}
