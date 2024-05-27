use ark_ff::PrimeField;

pub trait MultilinearPolynomialTrait<F: PrimeField> {
    fn partial_eval(&self, x: &Vec<(usize, F)>) -> Self;
    fn evaluate(&self, x: &Vec<(usize, F)>) -> F;
    fn number_of_vars(&self) -> usize;
    fn to_bytes(&self) -> Vec<u8>;
    fn relabel(&self) -> Self;
    fn additive_identity() -> Self;
    fn sum_over_the_boolean_hypercube(&self) -> F;
}