use std::marker::PhantomData;

use ark_ff::PrimeField;
use fiat_shamir_transcript::Transcript;
use polynomials::multilinear_polynomial::{eval_form::MLE, traits::MultilinearPolynomialTrait};

#[derive(Debug)]
pub struct SumcheckProof<F: PrimeField> {
    pub init_poly: MLE<F>,
    pub round_polys: Vec<MLE<F>>,
    pub claimed_sum: F
}

impl <F: PrimeField>SumcheckProof<F> {
    pub fn new(init_poly: MLE<F>, round_polys: Vec<MLE<F>>, claimed_sum: F) -> Self {
        Self { init_poly, round_polys, claimed_sum }
    }
}

pub struct Prover<F: PrimeField> {
    pub _marker: PhantomData<F>
}

impl <F: PrimeField>Prover<F> {

    pub fn prove(init_poly: MLE<F>) -> SumcheckProof<F> {

        let claimed_sum = init_poly.sum_over_the_boolean_hypercube();

        let mut transcript = Transcript::new();
        transcript.append(&init_poly.to_bytes());

        let mut round_polys: Vec<MLE<F>> = vec![];

        let mut initial_poly = init_poly.clone();
        
        for _ in 0..init_poly.num_of_vars {

            let round_poly = initial_poly.skip_one_and_sum_over_the_boolean_hypercube();

            let challenge = transcript.sample_field_element();
            
            initial_poly = initial_poly.partial_eval(&vec![(1, challenge)]);
            
            transcript.append(&round_poly.to_bytes());

            round_polys.push(round_poly);

        }

        SumcheckProof::new(init_poly, round_polys, claimed_sum)
    }
}


#[cfg(test)]
pub mod test {
    use polynomials::multilinear_polynomial::eval_form::MLE;
    use ark_bn254::Fq;

    use super::Prover;


    #[test]
    pub fn test_prove_sumcheck() {

        let val = vec![
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(2),
            Fq::from(2),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(2),
            Fq::from(2),
        ];

        let init_poly = MLE::new(&val);

        let proof = Prover::prove(init_poly);

        assert!(proof.claimed_sum == Fq::from(8), "Incorrect sum over the boolean hypercube");

        assert!(proof.round_polys[0].val == vec![
            Fq::from(4),
            Fq::from(4),
        ], "Incorrect round 1 poly");


    }
}
