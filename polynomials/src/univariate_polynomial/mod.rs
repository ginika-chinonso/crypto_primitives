use ark_ff::{BigInteger, PrimeField};
use ark_serialize::*;
use std::{
    collections::HashMap,
    fmt::Display,
    ops::{Add, Div, Mul},
};

// Univariate Polynomial
#[derive(Debug, Clone, PartialEq, CanonicalSerialize, CanonicalDeserialize)]
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
    pub fn interpolate(x_values: &Vec<F>, y_values: &Vec<F>) -> Self {
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
        res_poly.truncate()
    }

    // Checks if polynomial is a zero polynomial
    pub fn is_zero(&self) -> bool {
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

    pub fn additive_identity() -> Self {
        Self {
            coefficients: vec![F::zero()],
        }
    }

    pub fn multiplicative_identity() -> Self {
        Self {
            coefficients: vec![F::one()],
        }
    }

    pub fn to_bytes(&self) -> Vec<u8> {
        let mut res = Vec::new();
        for coefficient in self.coefficients.to_vec() {
            res.extend(coefficient.into_bigint().to_bytes_be())
        }
        res
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

// Implement native addition for univariate polynomial
impl<F: PrimeField> Add for UnivariatePolynomial<F> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self {
        if self.is_zero() {
            return rhs.clone();
        }

        if rhs.is_zero() {
            return self.clone();
        }

        let (mut longer, shorter) = if self.coefficients.len() >= rhs.coefficients.len() {
            (self.coefficients.clone(), &rhs.coefficients)
        } else {
            (rhs.coefficients.clone(), &self.coefficients)
        };

        for i in 0..shorter.len() {
            longer[i] += shorter[i];
        }

        Self::new(longer)
    }
}

impl<F: PrimeField> Display for UnivariatePolynomial<F> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for i in 0..self.coefficients.len() {
            if self.coefficients[i] == F::zero() {
                continue;
            }
            if i == 0 {
                write!(f, "{}", self.coefficients[i])?
            } else if self.coefficients[i] == F::one() {
                write!(f, " + x^{}", i)?
            } else {
                write!(f, " + {}x^{}", self.coefficients[i], i)?
            }
        }
        Ok(())
    }
}

// Implement native division for univariate polynomials
impl<F: PrimeField> Div for UnivariatePolynomial<F> {
    // Result should be the quotient and remainder

    type Output = (UnivariatePolynomial<F>, UnivariatePolynomial<F>);

    fn div(self, rhs: Self) -> Self::Output {
        assert!(
            self.coefficients.len() >= rhs.coefficients.len(),
            "Dividend degree should be higher than divisor"
        );
        let mut res = vec![];

        let mut remainder: Vec<F> = self.coefficients.into_iter().rev().collect();
        let divisor: Vec<F> = rhs.coefficients.into_iter().rev().collect();

        for i in 0..remainder.len() {
            if remainder.len() - i < divisor.len() {
                break;
            }

            let mut divisor_index = 0;
            let quotient = remainder[i] / divisor[divisor_index];

            for _ in 0..divisor.len() {
                remainder[i + divisor_index] =
                    remainder[i + divisor_index] - (quotient * divisor[divisor_index]);

                divisor_index += 1;
            }

            res.push(quotient);
        }

        res.reverse();
        remainder.reverse();

        // Returns the quotient and remainder polynomials
        (
            UnivariatePolynomial::new(res),
            UnivariatePolynomial::new(remainder),
        )
    }
}

// TODO:
// Implement FFT and IFFT

//////////////////////////////////////////////
/// TESTS
/// /////////////////////////////////////////
#[cfg(test)]
mod tests {

    use crate::multilinear_polynomial::{
        coef_form::{MultilinearMonomial, MultilinearPolynomial},
        traits::MultilinearPolynomialTrait,
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
    fn test_add_poly_buggy_test() {
        let coeffs1 = vec![Fq::from(1), Fq::from(16), Fq::from(13)];
        let coeffs2 = vec![Fq::from(16), Fq::from(9)];
        let poly1 = UnivariatePolynomial::new(coeffs1);
        let poly2 = UnivariatePolynomial::new(coeffs2);
        let poly_res = poly1 + poly2;
        let res_coeffs = vec![Fq::from(17), Fq::from(25), Fq::from(13)];
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
        let poly_res = UnivariatePolynomial::interpolate(&x_vals, &y_vals).truncate();
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
        let new_poly = multilinear_poly.to_univariate().unwrap();
        assert!(new_poly == UnivariatePolynomial::new(vec![Fq::from(14), Fq::from(9)]));
    }

    #[test]
    fn test_display_univariate_poly() {
        let poly = UnivariatePolynomial::new(vec![
            Fq::from(0),
            Fq::from(3),
            Fq::from(1),
            Fq::from(4),
            Fq::from(1),
            Fq::from(6),
        ]);

        println!("{}", poly);
    }

    #[test]
    pub fn test_poly_division() {
        let dividend = UnivariatePolynomial::new(vec![Fq::from(6), Fq::from(5), Fq::from(1)]);
        let divisor = UnivariatePolynomial::new(vec![Fq::from(2), Fq::from(1)]);

        let (quotient, remainder): (UnivariatePolynomial<Fq>, UnivariatePolynomial<Fq>) =
            dividend / divisor;

        assert!(remainder.truncate().is_zero(), "No remainder expected");
        assert!(
            quotient == UnivariatePolynomial::new(vec![Fq::from(3), Fq::from(1)]),
            "Incorrect quotient from division"
        );

        /////////////////

        let poly1 = UnivariatePolynomial::new(vec![Fq::from(15), Fq::from(8), Fq::from(1)]);
        let poly2 = UnivariatePolynomial::new(vec![Fq::from(5), Fq::from(1)]);
        let expected = UnivariatePolynomial::new(vec![Fq::from(3), Fq::from(1)]);

        let (res1quo, res1rem) = poly1 / poly2;
        assert!(res1rem.truncate().is_zero(), "No remainder expected");
        assert!(res1quo == expected, "Incorrect quotient from division");
    }
}
