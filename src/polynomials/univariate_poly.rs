use ark_ff::PrimeField;
use serde::{Deserialize, Serialize};
use std::{
    collections::HashMap,
    ops::{Add, Mul},
};

use super::multilinear_poly::MultilinearPolynomial;

// Univariate Polynomial
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct UnivariatePolynomial<F: PrimeField> {
    // Where elements index is their x power
    pub coefficients: Vec<F>,
}

// Polynomial implementation
impl<F: PrimeField> UnivariatePolynomial<F> {
    // Create new univariate polynomial
    pub fn new(coefficients: Vec<F>) -> Self {
        Self { coefficients }
    }

    // Evaluate polynomial at a point
    pub fn evaluate(&self, x: F) -> F {
        self.coefficients
            .iter()
            .rev()
            .fold(F::zero(), |res, val| (res * x) + val)
    }

    // Interpolate polynomial given x and y values
    pub fn interpolate(x_values: Vec<F>, y_values: Vec<F>) -> Self {
        assert_eq!(x_values.len(), y_values.len());
        let mut res_poly = Self::new(vec![]);
        let mut cache = HashMap::<(usize, usize), F>::new();
        for i in 0..x_values.len() {
            let mut y_poly = Self::new(vec![y_values[i]]);
            for j in 0..x_values.len() {
                if x_values[i] == x_values[j] {
                    continue;
                };
                let num = Self::new(vec![-x_values[j], F::one()]);

                let denom = if i > j {
                    Self::new(vec![cache.get(&(j, i)).unwrap().neg()])
                } else {
                    let r = (x_values[i] - x_values[j]).inverse().unwrap();
                    cache.insert((i, j), r.clone());
                    Self::new(vec![r])
                };
                // without caching
                // match (x_values[i] - x_values[j]).inverse() {
                //     Option::Some(val) => {
                //         denom = Polynomial::new(vec![val]);
                //     }
                //     Option::None => {
                //         panic!("Invalid inverse");
                //     }
                // }
                y_poly = y_poly * num * denom;
            }
            res_poly = res_poly + y_poly;
        }
        res_poly
    }

    // Checks if polynomial is a zero polynomial
    pub fn is_zero(self) -> bool {
        if self.coefficients.len() == 0 {
            true
        } else {
            false
        }
    }

    // removes zero terms from a polynomial
    pub fn truncate(mut self) -> Self {
        match self.coefficients.pop() {
            Option::Some(val) => {
                if val == F::zero() {
                    self.truncate()
                } else {
                    self.coefficients.push(val);
                    self
                }
            }
            Option::None => self,
        }
    }
}

// Implement native multiplication for univariate polynomial
impl<F: PrimeField> Mul for UnivariatePolynomial<F> {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        let mut res_array = vec![F::zero(); self.coefficients.len() + rhs.coefficients.len() - 1];

        let mut i = 0;

        loop {
            if i >= self.coefficients.len() {
                break;
            }
            let mut j = 0;
            loop {
                if j >= rhs.coefficients.len() {
                    break;
                };
                let val = self.coefficients[i] * rhs.coefficients[j];
                res_array[i + j] += val;
                j += 1;
            }
            i += 1;
        }

        UnivariatePolynomial::new(res_array)
    }
}

// Converts a multilinear polynomial to a univariate polynomial
impl<F: PrimeField> TryFrom<MultilinearPolynomial<F>> for UnivariatePolynomial<F> {
    type Error = &'static str;

    fn try_from(value: MultilinearPolynomial<F>) -> Result<Self, Self::Error> {
        let mut res = UnivariatePolynomial::<F>::new(vec![]);

        if value.terms.len() == 0 {
            res = UnivariatePolynomial {
                coefficients: vec![],
            };
        } else if value.terms[0].vars.len() > 1 {
            return Err("Not a univariate poly, try relabelling");
        } else if value.terms[0].vars.is_empty() {
            res = UnivariatePolynomial::<F>::new(vec![value.terms[0].coefficient]);
        } else {
            if value.terms[0].vars[0] == false {
                res = UnivariatePolynomial::<F>::new(vec![
                    value.terms[0].coefficient,
                    value
                        .terms
                        .get(1)
                        .map(|a| a.coefficient)
                        .unwrap_or(F::zero()),
                ]);
            } else {
                res = UnivariatePolynomial::<F>::new(vec![
                    value
                        .terms
                        .get(1)
                        .map(|a| a.coefficient)
                        .unwrap_or(F::zero()),
                    value.terms[0].coefficient,
                ]);
            }
        }

        Ok(res)
    }
}

// Implement native addition for univariate polynomial
impl<F: PrimeField> Add for UnivariatePolynomial<F> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self {
        let mut res_array = vec![];

        let mut i = 0;

        if self.coefficients.len() <= rhs.coefficients.len() {
            loop {
                if i == self.coefficients.len() {
                    break;
                }
                res_array.push(self.coefficients[i] + rhs.coefficients[i]);
                i += 1;
            }
            res_array.append(&mut rhs.coefficients[i..].to_vec());
        } else {
            loop {
                if i == rhs.coefficients.len() {
                    break;
                }
                res_array.append(&mut self.coefficients[i..].to_vec());
            }
        };

        Self::new(res_array)
    }
}

//////////////////////////////////////////////
/// TESTS
/// /////////////////////////////////////////
#[cfg(test)]
mod tests {
    use std::vec;

    use crate::polynomials::multilinear_poly::{
        MultilinearMonomial, MultilinearPolynomial,
    };

    use super::UnivariatePolynomial;
    use ark_ff::{Fp64, MontBackend, MontConfig};

    #[derive(MontConfig)]
    #[modulus = "17"]
    #[generator = "3"]
    pub struct FqConfig;
    pub type Fq = Fp64<MontBackend<FqConfig, 1>>;

    #[test]
    fn test_eval_poly() {
        let coeffs = vec![Fq::from(0), Fq::from(2)];
        let new_polynomial = UnivariatePolynomial::new(coeffs);
        let res = new_polynomial.evaluate(Fq::from(4));
        assert_eq!(res, Fq::from(8));
    }

    #[test]
    fn test_add_poly() {
        let coeffs1 = vec![Fq::from(1), Fq::from(2), Fq::from(3)];
        let coeffs2 = vec![Fq::from(2), Fq::from(3), Fq::from(4)];
        let poly1 = UnivariatePolynomial::new(coeffs1);
        let poly2 = UnivariatePolynomial::new(coeffs2);
        let poly_res = poly1 + poly2;
        let res_coeffs = vec![Fq::from(3), Fq::from(5), Fq::from(7)];
        assert_eq!(poly_res, UnivariatePolynomial::new(res_coeffs));
    }

    #[test]
    fn test_mul_poly() {
        let coeffs1 = vec![Fq::from(0), Fq::from(2), Fq::from(6)];
        let coeffs2 = vec![Fq::from(1), Fq::from(0), Fq::from(1)];
        let poly1 = UnivariatePolynomial::new(coeffs1);
        let poly2 = UnivariatePolynomial::new(coeffs2);
        let poly_res = poly1 * poly2;
        let res_coeffs = vec![
            Fq::from(0),
            Fq::from(2),
            Fq::from(6),
            Fq::from(2),
            Fq::from(6),
        ];
        assert_eq!(poly_res, UnivariatePolynomial::new(res_coeffs));
    }

    #[test]
    fn test_interpolate_poly() {
        let x_vals = vec![
            Fq::from(1),
            Fq::from(2),
            Fq::from(3),
            Fq::from(4),
            Fq::from(5),
        ];
        let y_vals = vec![
            Fq::from(6),
            Fq::from(11),
            Fq::from(18),
            Fq::from(27),
            Fq::from(38),
        ];
        let poly_res = UnivariatePolynomial::interpolate(x_vals, y_vals).truncate();
        let res_vals = vec![Fq::from(3), Fq::from(2), Fq::from(1)];
        assert_eq!(poly_res, UnivariatePolynomial::new(res_vals));
    }

    #[test]
    fn test_is_zero_poly() {
        let test_poly1: UnivariatePolynomial<Fq> = UnivariatePolynomial::new(vec![Fq::from(3)]);
        let test_poly2: UnivariatePolynomial<Fq> = UnivariatePolynomial::new(vec![]);
        assert_eq!(test_poly1.is_zero(), false);
        assert_eq!(test_poly2.is_zero(), true);
    }

    #[test]
    fn test_truncate_poly() {
        let test_poly = UnivariatePolynomial::new(vec![
            Fq::from(3),
            Fq::from(0),
            Fq::from(1),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
        ]);
        assert_eq!(
            test_poly.truncate(),
            UnivariatePolynomial::new(vec![Fq::from(3), Fq::from(0), Fq::from(1)])
        );
    }

    #[test]
    fn test_tryfrom() {
        let multilinear_poly = MultilinearPolynomial::new(vec![
            MultilinearMonomial::new(Fq::from(9), vec![true]),
            MultilinearMonomial::new(Fq::from(14), vec![false]),
        ]);
        let new_poly = UnivariatePolynomial::try_from(multilinear_poly).unwrap();
        assert!(new_poly == UnivariatePolynomial::new(vec![Fq::from(14), Fq::from(9)]));
    }
}
