use std::{error::Error, marker::PhantomData};

use ark_ff::PrimeField;
use fiat_shamir_transcript::Transcript;
use polynomials::{
    multilinear_polynomial::{
        traits::MultilinearPolynomialTrait, universal_mle::universal_mle::UniversalMLE,
    },
    univariate_polynomial::UnivariatePolynomial,
};

#[derive(Debug)]
pub struct SumcheckProof<F: PrimeField> {
    pub round_polys: Vec<UnivariatePolynomial<F>>,
    pub claimed_sum: F,
}

impl<F: PrimeField> SumcheckProof<F> {
    pub fn new(round_polys: Vec<UnivariatePolynomial<F>>, claimed_sum: F) -> Self {
        Self {
            round_polys,
            claimed_sum,
        }
    }
}

pub struct Sumcheck<F> {
    pub _field: PhantomData<F>,
}

impl<F: PrimeField> Sumcheck<F> {
    pub fn prove(initial_poly: UniversalMLE<F>) -> SumcheckProof<F> {
        let claimed_sum = initial_poly.sum_over_the_boolean_hypercube();

        let mut transcript = Transcript::new();
        // transcript.append(&initial_poly.to_bytes());

        let mut round_polys = vec![];

        let mut poly = initial_poly.clone();

        for _ in 0..initial_poly.number_of_vars() {
            let round_poly = poly.skip_one_and_sum_over_the_boolean_hypercube();

            transcript.append(&round_poly.to_bytes());
            let challenge: F = transcript.sample_field_element();

            poly = poly.partial_eval(&vec![(1, challenge)]);

            round_polys.push(round_poly);
        }

        SumcheckProof::new(round_polys, claimed_sum)
    }

    pub fn verify(
        initial_poly: UniversalMLE<F>,
        proof: SumcheckProof<F>,
    ) -> Result<bool, Box<dyn Error>> {
        let (challenges, _, claimed_sum) = Sumcheck::verify_partial(&proof)?;

        if initial_poly.evaluate(&challenges) != claimed_sum {
            return Err(format!("Verifier oracle query check failed").into());
        }

        Ok(true)
    }

    pub fn verify_partial(
        proof: &SumcheckProof<F>,
    ) -> Result<(Vec<(usize, F)>, Transcript, F), Box<dyn Error>> {
        let mut transcript = Transcript::new();

        let mut claimed_sum = proof.claimed_sum;
        let mut challenges = vec![];

        for i in 0..proof.round_polys.len() {
            if claimed_sum
                != proof.round_polys[i].evaluate(F::zero())
                    + proof.round_polys[i].evaluate(F::one())
            {
                return Err(format!("Verification for round {} failed", i).into());
            }

            transcript.append(&proof.round_polys[i].to_bytes());
            let challenge: F = transcript.sample_field_element();

            claimed_sum = proof.round_polys[i].evaluate(challenge);

            challenges.push((i + 1, challenge));
        }

        Ok((challenges, transcript, claimed_sum))
    }
}

#[cfg(test)]
pub mod test {
    use ark_bn254::Fq;
    use polynomials::multilinear_polynomial::{
        eval_form::MLE,
        universal_mle::{ops::Ops, universal_mle::UniversalMLE},
    };

    use super::Sumcheck;

    pub fn create_universal_mle() -> UniversalMLE<Fq> {
        // 2bc * 3bc * 5bc = 30b^3c^3
        let poly_1 = UniversalMLE::Value(MLE::new(&vec![
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(2),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(2),
        ]));
        let poly_2 = UniversalMLE::Value(MLE::new(&vec![
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(3),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(3),
        ]));
        let poly_3 = UniversalMLE::Value(MLE::new(&vec![
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(5),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(5),
        ]));

        UniversalMLE::new(vec![poly_1, poly_2, poly_3], Ops::new().mul)
    }

    #[test]
    pub fn test_universal_mle_sumcheck() {
        let poly = create_universal_mle();

        let proof = Sumcheck::prove(poly.clone());

        let verify = Sumcheck::verify(poly, proof);

        assert!(
            verify.expect("Sumcheck verification failed"),
            "Sumcheck verification falied"
        );
    }
}
