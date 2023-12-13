use ark_ff::{BigInteger, PrimeField};
use serde::{Deserialize, Serialize};

use crate::{
    fiat_shamir_transcript::transcript::Transcript,
    polynomials::{multilinear_poly::MultilinearPolynomial, univariate_poly::UnivariatePolynomial},
};

use super::prover::Prover;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Proof<F: PrimeField> {
    initial_poly: MultilinearPolynomial<F>,
    sum: F,
    rounds_poly: Vec<UnivariatePolynomial<F>>,
}

impl<F: PrimeField> Proof<F> {
    pub fn new(initial_poly: MultilinearPolynomial<F>, sum: F) -> Self {
        Self {
            initial_poly,
            sum,
            rounds_poly: vec![],
        }
    }
}

pub struct Sumcheck {}

impl Sumcheck {
    pub fn prove<F: PrimeField>(poly: MultilinearPolynomial<F>, sum: F) -> Proof<F> {
        let mut proof = Proof::new(poly.clone(), sum);
        let mut transcript = Transcript::new();
        let mut challenges: Vec<F> = vec![];

        transcript.add_multivariate_poly(&poly);
        transcript.append(sum.into_bigint().to_bytes_be().as_slice());

        let prover = Prover::<F>::new(poly.clone(), sum);

        for _ in 0..poly.terms[0].vars.len() {
            // if i == 0 {
            //     let round_poly = prover.prove(&[]);
            //     proof.rounds_poly.push(round_poly);
            // };
            let round_poly = prover.prove(&challenges);
            transcript.add_univariate_poly(&round_poly);
            proof.rounds_poly.push(round_poly);

            challenges.push(transcript.sample_field_element());
        }

        proof
    }

    pub fn verify<F: PrimeField>(proof: Proof<F>) -> bool {
        let mut transcript = Transcript::new();
        let mut challenges = vec![];

        transcript.add_multivariate_poly(&proof.initial_poly);
        transcript.append(proof.sum.into_bigint().to_bytes_be().as_slice());

        for i in 0..proof.initial_poly.terms[0].vars.len() {
            let verifier_check =
                proof.rounds_poly[i].evaluate(F::zero()) + proof.rounds_poly[i].evaluate(F::one());

            if i == 0 {
                if verifier_check != proof.sum {
                    return false;
                };
            };

            transcript.add_univariate_poly(&proof.rounds_poly[i]);
            challenges.push(transcript.sample_field_element::<F>());

            let last_round_sum = if i == 0 {
                proof.sum
            } else {
                proof.rounds_poly[i - 1].evaluate(challenges[i])
            };

            if challenges.len() == proof.initial_poly.clone().terms[0].vars.len() {
                return last_round_sum
                    == proof.initial_poly.evaluate(
                        challenges
                            .iter()
                            .enumerate()
                            .map(|(a, b)| (a, *b))
                            .collect(),
                    );
            }
        }

        true
    }
}

#[cfg(test)]
mod tests {
    use crate::polynomials::multilinear_poly::{MultilinearMonomial, MultilinearPolynomial};

    use super::Sumcheck;

    use ark_ff::{Fp64, MontBackend, MontConfig};

    #[derive(MontConfig)]
    #[modulus = "17"]
    #[generator = "3"]
    pub struct FqConfig;
    pub type Fq = Fp64<MontBackend<FqConfig, 1>>;

    #[test]
    fn test_noninteractive_sumcheck() {
        let proof = Sumcheck::prove(
            MultilinearPolynomial::new(vec![
                MultilinearMonomial::new(Fq::from(2), vec![true, true, false]),
                MultilinearMonomial::new(Fq::from(3), vec![false, true, true]),
            ]),
            Fq::from(10),
        );

        dbg!(&proof);
        let accepted = Sumcheck::verify(proof);

        assert!(accepted);
    }
}
