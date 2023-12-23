use ark_ff::{BigInteger, PrimeField};
use serde::{Deserialize, Serialize};

use crate::{
    fiat_shamir_transcript::transcript::Transcript,
    polynomials::multilinear_poly::MultilinearPolynomialTrait,
};

use super::prover::Prover;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Proof<F: PrimeField, MPT: MultilinearPolynomialTrait<F> + Clone> {
    initial_poly: MPT,
    sum: F,
    rounds_poly: Vec<MPT>,
}

impl<F: PrimeField, MPT: MultilinearPolynomialTrait<F> + Clone> Proof<F, MPT> {
    pub fn new(initial_poly: MPT, sum: F) -> Self {
        Self {
            initial_poly,
            sum,
            rounds_poly: vec![],
        }
    }
}

pub struct Sumcheck {}

impl Sumcheck {
    pub fn prove<
        F: PrimeField,
        MPT: MultilinearPolynomialTrait<F> + Clone + std::ops::Add<Output = MPT>,
    >(
        poly: MPT,
        sum: F,
    ) -> Proof<F, MPT> {
        let mut proof = Proof::new(poly.clone(), sum);
        let mut transcript = Transcript::new();
        let mut challenges: Vec<F> = vec![];

        transcript.add_multivariate_poly(&poly);
        transcript.append(sum.into_bigint().to_bytes_be().as_slice());

        let prover = Prover::<F, MPT>::new(poly.clone(), sum);

        for _ in 0..poly.number_of_vars() {
            // if i == 0 {
            //     let round_poly = prover.prove(&[]);
            //     proof.rounds_poly.push(round_poly);
            // };
            let round_poly = prover.prove(&challenges);
            transcript.add_multivariate_poly(&round_poly);
            proof.rounds_poly.push(round_poly);

            challenges.push(transcript.sample_field_element());
            dbg!(&challenges);
        }

        proof
    }

    pub fn verify<F: PrimeField, MPT: MultilinearPolynomialTrait<F> + Clone + std::fmt::Debug>(
        proof: Proof<F, MPT>,
    ) -> bool {
        let mut transcript = Transcript::new();
        let mut challenges: Vec<(usize, F)> = vec![];
        transcript.add_multivariate_poly(&proof.initial_poly);
        transcript.append(proof.sum.into_bigint().to_bytes_be().as_slice());

        for i in 0..proof.initial_poly.number_of_vars() {
            let verifier_check = proof.rounds_poly[i].evaluate(vec![(0, F::zero())])
                + proof.rounds_poly[i].evaluate(vec![(0, F::one())]);

            if i == 0 {
                if verifier_check != proof.sum {
                    return false;
                };
            };

            transcript.add_multivariate_poly(&proof.rounds_poly[i]);
            challenges.extend(vec![(i, transcript.sample_field_element::<F>())]);

            let last_round_sum = if i == 0 {
                proof.sum
            } else {
                let (_, challenge) = challenges[i];
                proof.rounds_poly[i].evaluate(vec![(0, challenge)])
            };

            if challenges.len() == proof.initial_poly.number_of_vars() {
                let res = proof
                    .initial_poly
                    .evaluate(challenges.into_iter().collect());

                return last_round_sum == res;
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
        let accepted = Sumcheck::verify(proof);

        assert!(accepted);
    }
}
