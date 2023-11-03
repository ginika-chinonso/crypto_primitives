use ark_ff::Field;
use std::ops::{Add, Mul};

#[derive(Debug, Clone, PartialEq)]
pub struct Polynomial<F: Field> {
    // Where elements index is their x power
    pub coefficients: Vec<F>,
}

pub trait PolynomialTrait<F: Field> {
    fn new(coefficients: Vec<F>) -> Polynomial<F>;
    fn evaluate(&self, x: F) -> F;
    fn interpolate(x_values: Vec<F>, y_values: Vec<F>) -> Polynomial<F>;
    fn is_zero(self) -> bool;
    fn truncate(self) -> Polynomial<F>;
    // fn sparse_to_dense(self) -> Polynomial<F>;
}

impl<F: Field + std::convert::From<i32>> PolynomialTrait<F> for Polynomial<F> {
    fn new(coefficients: Vec<F>) -> Polynomial<F> {
        Polynomial { coefficients }
    }

    fn evaluate(&self, x: F) -> F {
        self.coefficients
            .iter()
            .rev()
            .fold(F::zero(), |res, val| (res * x) + val)
    }

    fn interpolate(x_values: Vec<F>, y_values: Vec<F>) -> Polynomial<F> {
        assert_eq!(x_values.len(), y_values.len());
        let mut res_poly = Polynomial::new(vec![]);
        for i in 0..x_values.len() {
            let mut y_poly = Polynomial::new(vec![y_values[i]]);
            for j in 0..x_values.len() {
                if x_values[i] == x_values[j] {
                    continue;
                };
                let num = Polynomial::new(vec![-x_values[j], F::from(1)]);
                let mut denom = Polynomial::new(vec![F::one()]);

                // match (x_values[i] - x_values[j]).inverse() {
                //     Option::Some(val) => {
                //         denom = Polynomial::new(vec![val]);
                //     }
                //     Option::None => {
                //         panic!("Invalid inverse");
                //     }
                // }
                // y_poly = y_poly * num * denom;
                dbg!("First mul");
                // let n_poly = num * denom;
                // dbg!("Second mul");
                y_poly.clone() * num * denom;
            }
            res_poly = res_poly + y_poly;
        }
        // dbg!(val_store);
        res_poly
    }

    fn is_zero(self) -> bool {
        if self.coefficients.len() == 0 {
            true
        } else {false}
    }

    fn truncate(mut self) -> Polynomial<F> {
        match self.coefficients.pop() {
            Option::Some(val) => {
                if val == F::zero() {
                    self.truncate()
                } else {
                    self.coefficients.push(val);
                    self
                }
            },
            Option::None => {self}
        }
    }
}

impl<F: Field + std::convert::From<i32>> Mul for Polynomial<F> {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        dbg!("Creating res array");
        let degree = self.coefficients.len() + rhs.coefficients.len() - 1;
        let mut res_array = vec![F::zero(); degree];
        // let mut res_array = vec![F::zero()];
        dbg!("Done creating");
        dbg!(self.coefficients.len() + rhs.coefficients.len() - 1);

        // let mut i = 0;

        // loop {
        //     if i >= self.coefficients.len() {
        //         break;
        //     }
        //     let mut j = 0;
        //     loop {
        //         if j >= rhs.coefficients.len() {
        //             break;
        //         };
        //         let val = self.coefficients[i] * rhs.coefficients[j];
        //         res_array[i + j] += val;
        //         j += 1;
        //     }
        //     i += 1;
        // }
            dbg!("Creating new polynomial");
        let p = Polynomial::new(res_array);
        dbg!("Poly created");
        p
    }
}

impl<F: Field + std::convert::From<i32>> Add for Polynomial<F> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self {
        let mut res_array = vec![];

        // loop {
        //     if &self.coefficients.len() == &rhs.coefficients.len() {
        //         break;
        //     }
        //     if &self.coefficients.len() < &rhs.coefficients.len() {
        //         self.coefficients.push(F::zero());
        //     } else {
        //         rhs.coefficients.push(F::zero());
        //     }
        // }

        // let mut i = 0;

        // loop {
        //     if i == self.coefficients.len(){
        //         break;
        //     }
        //     res_array.push(self.coefficients[i] + rhs.coefficients[i]);
        //     i += 1;
        // }

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

#[cfg(test)]
mod tests {
    use std::vec;

    use super::{Polynomial, PolynomialTrait};
    use ark_ff::{Fp64, MontBackend, MontConfig};

    #[derive(MontConfig)]
    #[modulus = "17"]
    #[generator = "3"]
    pub struct FqConfig;
    pub type Fq = Fp64<MontBackend<FqConfig, 1>>;

    #[test]
    fn test_eval_poly() {
        let coeffs = vec![Fq::from(0), Fq::from(2)];
        let new_polynomial = Polynomial::new(coeffs);
        let res = new_polynomial.evaluate(Fq::from(4));
        assert_eq!(res, Fq::from(8));
    }

    #[test]
    fn test_add_poly() {
        let coeffs1 = vec![Fq::from(1), Fq::from(2), Fq::from(3)];
        let coeffs2 = vec![Fq::from(2), Fq::from(3), Fq::from(4)];
        let poly1 = Polynomial::new(coeffs1);
        let poly2 = Polynomial::new(coeffs2);
        let poly_res = poly1 + poly2;
        let res_coeffs = vec![Fq::from(3), Fq::from(5), Fq::from(7)];
        assert_eq!(poly_res, Polynomial::new(res_coeffs));
    }

    #[test]
    fn test_mul_poly() {
        let coeffs1 = vec![Fq::from(0), Fq::from(2), Fq::from(6)];
        let coeffs2 = vec![Fq::from(1), Fq::from(0), Fq::from(1)];
        let poly1 = Polynomial::new(coeffs1);
        let poly2 = Polynomial::new(coeffs2);
        let poly_res = poly1 * poly2;
        let res_coeffs = vec![
            Fq::from(0),
            Fq::from(2),
            Fq::from(6),
            Fq::from(2),
            Fq::from(6),
        ];
        assert_eq!(poly_res, Polynomial::new(res_coeffs));
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
        let poly_res = Polynomial::interpolate(x_vals, y_vals).truncate();
        let res_vals = vec![Fq::from(3), Fq::from(2), Fq::from(1)];
        assert_eq!(poly_res, Polynomial::new(res_vals));
    }

    #[test]
    fn test_is_zero_poly() {
        let test_poly1: Polynomial<Fq> = Polynomial::new(vec![Fq::from(3)]);
        let test_poly2: Polynomial<Fq> = Polynomial::new(vec![]);
        assert_eq!(test_poly1.is_zero(), false);
        assert_eq!(test_poly2.is_zero(), true);
    }

    #[test]
    fn test_truncate_poly() {
        let test_poly = Polynomial::new(vec![Fq::from(3), Fq::from(0), Fq::from(1), Fq::from(0), Fq::from(0), Fq::from(0)]);
        assert_eq!(test_poly.truncate(), Polynomial::new(vec![Fq::from(3), Fq::from(0), Fq::from(1)]));
    }
}

// benchmarking
