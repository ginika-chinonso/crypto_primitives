use std::marker::PhantomData;

use ark_ff::PrimeField;
use fiat_shamir_transcript::Transcript;
use polynomials::{
    multilinear_polynomial::{eval_form::MLE, traits::MultilinearPolynomialTrait},
    univariate_polynomial::UnivariatePolynomial,
};

use crate::universal_mle::universal_mle::SumcheckProof;

pub struct Prover<F: PrimeField> {
    pub _marker: PhantomData<F>,
}

impl<F: PrimeField> Prover<F> {
    pub fn prove_sumcheck(init_poly: MLE<F>) -> SumcheckProof<F> {
        let claimed_sum = init_poly.sum_over_the_boolean_hypercube();

        let mut transcript = Transcript::new();
        transcript.append(&init_poly.to_bytes());

        let mut round_polys: Vec<UnivariatePolynomial<F>> = vec![];

        let mut initial_poly = init_poly.clone();

        for _ in 0..init_poly.num_of_vars {
            let round_poly = initial_poly.skip_one_and_sum_over_the_boolean_hypercube(2);

            let challenge = transcript.sample_field_element();

            initial_poly = initial_poly.partial_eval(&vec![(1, challenge)]);

            transcript.append(&round_poly.to_bytes());

            round_polys.push(round_poly);
        }

        SumcheckProof::new(round_polys, claimed_sum)
    }

    pub fn prove_sum_of_product(f_poly: &MLE<F>, g_poly: &MLE<F>) -> SumcheckProof<F> {
        let mut transcript = Transcript::new();
        transcript.append(&f_poly.to_bytes());
        transcript.append(&g_poly.to_bytes());

        let mut round_polys = vec![];

        assert_eq!(
            f_poly.val.len(),
            g_poly.val.len(),
            "Both f and g polynomials should have the same length"
        );

        let claimed_sum = f_poly.sum_product_over_the_boolean_hypercube(&vec![g_poly.clone()]);

        let mut f = f_poly.clone();
        let mut g = g_poly.clone();

        for _ in 0..f_poly.num_of_vars {
            let round_poly =
                f.skip_one_and_sum_product_over_the_boolean_hypercube(&vec![g.clone()]);

            let challenge = transcript.sample_field_element();

            transcript.append(&round_poly.to_bytes());

            round_polys.push(round_poly);

            f = f.partial_eval(&vec![(1, challenge)]);
            g = g.partial_eval(&vec![(1, challenge)]);
        }

        SumcheckProof::new(round_polys, claimed_sum)
    }
}

#[cfg(test)]
pub mod test {
    use ark_bn254::Fq;
    use polynomials::{
        multilinear_polynomial::eval_form::MLE, univariate_polynomial::UnivariatePolynomial,
    };

    use crate::eval_form::verifier::Verifier;

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

        let proof = Prover::prove_sumcheck(init_poly);

        assert!(
            proof.claimed_sum == Fq::from(8),
            "Incorrect sum over the boolean hypercube"
        );

        assert!(
            proof.round_polys[0].coefficients == vec![Fq::from(4)],
            "Incorrect round 1 poly"
        );
    }

    #[test]
    pub fn prove_sumcheck_product() {
        let val1 = vec![
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(2),
        ];

        let f_poly = MLE::new(&val1);

        let val2 = vec![
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(3),
        ];

        let g_poly = MLE::new(&val2);

        let proof = Prover::prove_sum_of_product(&f_poly, &g_poly);

        let verify = Verifier::verify_sumcheck_product(vec![f_poly, g_poly], proof);

        assert!(verify.unwrap(), "Invalid sumcheck proof");
    }

    #[test]
    pub fn test_interpolate() {
        let y_vals = vec![Fq::from(0), Fq::from(6), Fq::from(36)];
        let poly = UnivariatePolynomial::<Fq>::interpolate(
            &vec![Fq::from(0), Fq::from(1), Fq::from(2)],
            &y_vals,
        );
        dbg!(poly);
    }
}
