use ark_ff::PrimeField;
use polynomials::{
    multilinear_polynomial::traits::MultilinearPolynomialTrait,
    univariate_polynomial::UnivariatePolynomial,
};
use utils::get_binary_string;

#[derive(Debug, Clone)]
pub struct Transcript<F: PrimeField> {
    pub poly: UnivariatePolynomial<F>,
    pub random_challenge: F,
    pub sum: F,
}

#[derive(Debug, Clone)]
pub struct Prover<F: PrimeField, MPT: MultilinearPolynomialTrait<F> + Clone> {
    pub initial_poly: MPT,
    pub initial_sum: F,
    pub rounds_transcript: Vec<Transcript<F>>,
}

impl<F: PrimeField, MPT: MultilinearPolynomialTrait<F> + Clone + std::ops::Add<Output = MPT>>
    Prover<F, MPT>
{
    pub fn new(poly: MPT, sum: F) -> Self {
        Prover {
            initial_poly: poly.clone(),
            initial_sum: sum,
            rounds_transcript: vec![],
        }
    }

    pub fn prove(&self, challenges: &[F]) -> UnivariatePolynomial<F> {
        let mut num_of_vars = self.initial_poly.number_of_vars();
        let mut round_poly = MPT::additive_identity();

        num_of_vars -= challenges.len();
        num_of_vars -= 1;

        for i in 0..2_usize.pow(num_of_vars as u32) {
            let eval_points = get_binary_string(i, num_of_vars);
            let mut eval_points_as_field: Vec<(usize, F)> = eval_points
                .chars()
                .enumerate()
                .map(|(a, b)| {
                    if b == '0' {
                        (a + challenges.len() + 1, F::zero())
                    } else {
                        (a + challenges.len() + 1, F::one())
                    }
                })
                .collect();

            for i in 0..challenges.len() {
                eval_points_as_field.push((i, challenges[i]));
            }
            let eval_res = self
                .initial_poly
                .clone()
                .partial_eval(&eval_points_as_field);
            round_poly = round_poly + eval_res;
        }
        let relabeled_round_poly = round_poly.relabel();
        relabeled_round_poly.to_univariate().unwrap()
    }
}
