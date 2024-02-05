use ark_ff::{BigInteger, PrimeField};
use serde::{Deserialize, Serialize};

use crate::{
    fiat_shamir_transcript::transcript::Transcript,
    polynomials::multilinear_poly::MultilinearPolynomialTrait,
};

use super::prover::Prover;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SumcheckProof<F: PrimeField, MPT: MultilinearPolynomialTrait<F> + Clone> {
    pub sum: F,
    pub number_of_vars: usize,
    pub challenges: Vec<F>,
    pub rounds_poly: Vec<MPT>,
}

impl<F: PrimeField, MPT: MultilinearPolynomialTrait<F> + Clone> SumcheckProof<F, MPT> {
    pub fn new(sum: F, number_of_vars: usize) -> Self {
        Self {
            sum,
            number_of_vars,
            challenges: vec![],
            rounds_poly: vec![],
        }
    }

    pub fn to_bytes(&self) -> Vec<u8> {
        let mut res = vec![];

        res.extend(self.number_of_vars.to_be_bytes());
        res.extend(self.sum.into_bigint().to_bytes_be());
        let res = self.challenges.iter().fold(res, |mut res, challenge| {
            res.extend(challenge.into_bigint().to_bytes_be());
            res
        });

        self.rounds_poly.iter().fold(res, |mut res, poly| {
            res.extend(poly.to_bytes());
            res
        })
    }
}

pub struct Sumcheck {}

impl Sumcheck {
    pub fn prove<
        F: PrimeField,
        MPT: MultilinearPolynomialTrait<F> + Clone + std::ops::Add<Output = MPT> + std::fmt::Debug,
    >(
        initial_poly: &MPT,
        sum: &F,
    ) -> SumcheckProof<F, MPT> {
        let mut proof = SumcheckProof::new(*sum, initial_poly.number_of_vars());
        let mut transcript = Transcript::new();
        let mut challenges: Vec<F> = vec![];

        // transcript.add_multivariate_poly(&poly);
        transcript.append(sum.into_bigint().to_bytes_be().as_slice());

        let prover = Prover::<F, MPT>::new(initial_poly.clone(), *sum);

        for _ in 0..initial_poly.number_of_vars() {
            let round_poly = prover.prove(&challenges);
            transcript.add_multivariate_poly(&round_poly);
            proof.rounds_poly.push(round_poly);

            challenges.push(transcript.sample_field_element());
        }

        proof.challenges.extend(challenges);

        proof
    }

    pub fn verify<F: PrimeField, MPT: MultilinearPolynomialTrait<F> + Clone + std::fmt::Debug>(
        mut proof: SumcheckProof<F, MPT>,
        initial_poly: MPT,
    ) -> Result<bool, String> {
        let (challenges, _) = Sumcheck::verify_partial(&mut proof)?;

        let res = initial_poly.evaluate(&challenges);

        Ok(proof.sum == res)
    }

    // pub fn prove_partial<
    //     F: PrimeField,
    //     MPT: MultilinearPolynomialTrait<F> + Clone + std::ops::Add<Output = MPT>,
    // >(
    //     initial_poly: MPT,
    //     poly: MPT,
    //     sum: F,
    // ) -> Proof<F, MPT> {
    //     let mut proof = Proof::new(sum);
    //     let mut transcript = Transcript::new();
    //     let mut challenges: Vec<F> = vec![];

    //     transcript.add_multivariate_poly(&poly);
    //     transcript.append(sum.into_bigint().to_bytes_be().as_slice());

    //     let prover = Prover::<F, MPT>::new(poly.clone(), sum);

    //     for _ in 0..poly.number_of_vars() - 1{

    //         let round_poly = prover.prove(&challenges);
    //         transcript.add_multivariate_poly(&round_poly);
    //         proof.rounds_poly.push(round_poly);

    //         challenges.push(transcript.sample_field_element());
    //         dbg!(&challenges);
    //     }

    //     proof
    // }

    pub fn verify_partial<
        F: PrimeField,
        MPT: MultilinearPolynomialTrait<F> + Clone + std::fmt::Debug,
    >(
        proof: &mut SumcheckProof<F, MPT>,
    ) -> Result<(Vec<(usize, F)>, F), String> {
        let mut transcript = Transcript::new();
        let mut challenges: Vec<(usize, F)> = vec![];
        // transcript.add_multivariate_poly(&proof.initial_poly);
        transcript.append(proof.sum.into_bigint().to_bytes_be().as_slice());

        for i in 0..proof.number_of_vars {

            let verifier_check = proof.rounds_poly[i].evaluate(&vec![(0, F::zero())])
                + proof.rounds_poly[i].evaluate(&vec![(0, F::one())]);

            if verifier_check != proof.sum {
                return Err("Invalid proof".to_string());
            };

            transcript.add_multivariate_poly(&proof.rounds_poly[i]);
            challenges.extend(vec![(i, transcript.sample_field_element::<F>())]);

            let (_, challenge) = challenges[i];

            proof.sum = proof.rounds_poly[i].evaluate(&vec![(0, challenge)]);
        }
        // dbg!(&challenges);
        Ok((challenges, proof.sum))
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
            &MultilinearPolynomial::new(vec![
                MultilinearMonomial::new(Fq::from(2), vec![true, true, false]),
                MultilinearMonomial::new(Fq::from(3), vec![false, true, true]),
            ]),
            &Fq::from(10),
        );
        let accepted = Sumcheck::verify(
            proof,
            MultilinearPolynomial::new(vec![
                MultilinearMonomial::new(Fq::from(2), vec![true, true, false]),
                MultilinearMonomial::new(Fq::from(3), vec![false, true, true]),
            ]),
        );

        assert!(accepted.unwrap());
    }
}
