use crate::polynomials::{
    multilinear_poly::MultilinearPolynomial, univariate_poly::UnivariatePolynomial,
};
use ark_ff::PrimeField;

#[derive(Debug, Clone)]
pub struct Transcript<F: PrimeField> {
    pub poly: UnivariatePolynomial<F>,
    pub random_challenge: F,
    pub sum: F,
}

#[derive(Debug, Clone)]
pub struct Prover<F: PrimeField> {
    pub initial_poly: MultilinearPolynomial<F>,
    pub initial_sum: F,
    pub rounds_transcript: Vec<Transcript<F>>,
}

impl<F: PrimeField> Prover<F> {
    pub fn new(poly: MultilinearPolynomial<F>, sum: F) -> Self {
        Prover {
            initial_poly: poly.clone(),
            initial_sum: sum,
            rounds_transcript: vec![],
        }
    }

    pub fn prove(&self, challenges: &[F]) -> UnivariatePolynomial<F> {
        let mut num_of_vars = self.initial_poly.terms[0].vars.len();
        let mut round_poly = MultilinearPolynomial::<F>::new(vec![]);
        num_of_vars -= challenges.len();
        num_of_vars -= 1;

        for i in 0..2_usize.pow(num_of_vars as u32) {
            dbg!(&i);
            dbg!(&num_of_vars);
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
                .partial_eval(eval_points_as_field)
                .relabel();
            round_poly = round_poly.add(eval_res.clone());
        }
        dbg!(&round_poly);
        UnivariatePolynomial::try_from(round_poly.relabel()).unwrap()
    }
}

// Get padded binary string of a decimal number
fn get_binary_string(index: usize, max_bit_count: usize) -> String {
    if max_bit_count == 0 {
        return "".to_string();
    }
    let binary = format!("{:b}", index);
    "0".repeat(max_bit_count - binary.len()) + &binary
}
