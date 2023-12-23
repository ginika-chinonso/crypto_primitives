use crate::polynomials::multilinear_poly::{MultilinearMonomial, MultilinearPolynomial};
use ark_ff::PrimeField;

// Get padded binary string of a decimal number
pub fn get_binary_string(index: usize, max_bit_count: usize) -> String {
    if max_bit_count == 0 {
        return "".to_string();
    }
    let binary = format!("{:b}", index);
    "0".repeat(max_bit_count - binary.len()) + &binary
}

// Returns the multilinear polynomial representing a bit value
pub fn check_bit<F: PrimeField>(bit: usize) -> MultilinearPolynomial<F> {
    let mut res = MultilinearPolynomial::new(vec![]);
    // The check: xiwi + (1 - xi)(1 - wi)
    if bit == 0 {
        res = MultilinearPolynomial::new(vec![
            MultilinearMonomial::new(F::one(), vec![false]),
            MultilinearMonomial::new(F::one().neg(), vec![true]),
        ]);
    };
    if bit == 1 {
        res = MultilinearPolynomial::new(vec![MultilinearMonomial::new(F::one(), vec![true])]);
    }

    res
}
