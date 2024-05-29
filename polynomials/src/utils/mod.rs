use super::multilinear_polynomial::coef_form::{MultilinearMonomial, MultilinearPolynomial};
use ark_ff::PrimeField;

// Get padded binary string of a decimal number
pub fn get_binary_string(index: usize, max_bit_count: usize) -> String {
    if max_bit_count == 0 {
        return "".to_string();
    }
    let binary = format!("{:b}", index);
    "0".repeat(max_bit_count.checked_sub(binary.len()).unwrap_or(0)) + &binary
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

/// Maps a binary string to a number
pub fn selector_to_index(selector: &[bool]) -> usize {
    let mut sum = 0;
    let mut adder = 1;

    for i in 0..selector.len() {
        if selector[i] {
            sum += adder;
        }
        adder *= 2;
    }

    sum
}

pub fn get_sib(index: usize, num_of_vars: usize, var: usize) -> usize {
    index + (2_usize.pow(num_of_vars as u32) / 2_usize.pow(var as u32))
}

pub fn pad_to_len(num: usize, len: usize) -> String {
    let m = format!("{:b}", num);
    if m.len() < len {
        let n = String::from("0").repeat(len - m.len());
        return n + &m;
    };
    m
}

pub fn contains_var(val: String, var: usize) -> bool {
    if val.chars().nth(var).unwrap().to_string() == "0".to_string() {
        return false;
    } else {
        return true;
    }
}
