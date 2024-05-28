use crate::{
    multilinear_polynomial::traits::MultilinearPolynomialTrait,
    univariate_polynomial::UnivariatePolynomial, utils::get_sib,
};
use ark_ff::{BigInteger, PrimeField};
use ark_serialize::*;
use std::ops::{Add, Mul};

#[derive(Debug, Clone, CanonicalSerialize, CanonicalDeserialize)]
pub struct MLE<F: PrimeField> {
    pub num_of_vars: usize,
    // The val vector contains the evaluation of the mle over the boolean hypercube
    pub val: Vec<F>,
}

impl<F: PrimeField> MLE<F> {
    pub fn new(val: &Vec<F>) -> Self {
        // assert!(val.len().is_power_of_two(), "Number of evaluations should be a power of two");
        let num_of_vars = (val.len() as f64).log2().ceil() as usize;
        Self {
            num_of_vars,
            val: val.to_vec(),
        }
    }
}

impl<F: PrimeField> MultilinearPolynomialTrait<F> for MLE<F> {
    fn partial_eval(&self, points: &Vec<(usize, F)>) -> Self {
        let mut new_poly = self.clone();
        let mut res = vec![];
        let mut points = points.clone();
        points.sort();

        let mut num_of_vars = new_poly.num_of_vars;
        let mut previously_evaluated = 0;

        for i in 0..points.len() {
            let (var, val) = if previously_evaluated == 0 {
                previously_evaluated = points[i].0;
                points[i]
            } else if previously_evaluated < points[i].0 {
                previously_evaluated = points[i].0;
                (points[i].0 - i, points[i].1)
            } else {
                previously_evaluated = points[i].0;
                points[i]
            };

            assert!(var <= num_of_vars, "Variable not found");
            let mut index = 0;
            let mut shift = 2_usize.pow(num_of_vars as u32) / 2_usize.pow(var as u32);

            for _ in 0..(new_poly.val.len() / 2) {
                if shift == 0 {
                    shift = 2_usize.pow(num_of_vars as u32) / 2_usize.pow(var as u32);
                    index += shift;
                }

                let left = new_poly.val[index];

                let right = new_poly.val[get_sib(index, num_of_vars, var)];

                let new = (right * val) + ((F::ONE - val) * left);

                res.push(new);

                shift -= 1;
                index += 1;
            }

            new_poly = Self::new(&res);
            res = vec![];
            num_of_vars = new_poly.num_of_vars;
        }

        new_poly
    }

    fn evaluate(&self, points: &Vec<(usize, F)>) -> F {
        let res = self.partial_eval(&points);
        res.val[0]
    }

    fn sum_over_the_boolean_hypercube(&self) -> F {
        self.val.iter().fold(F::zero(), |acc, val| acc + val)
    }

    fn number_of_vars(&self) -> usize {
        self.num_of_vars
    }

    fn to_bytes(&self) -> Vec<u8> {
        let mut res = Vec::new();

        res.append(self.clone().num_of_vars.to_be_bytes().to_vec().as_mut());

        for value in &self.val {
            res.append(value.into_bigint().to_bytes_be().as_mut());
        }

        res
    }

    fn relabel(&self) -> Self {
        self.clone()
    }

    fn additive_identity() -> Self {
        Self {
            val: vec![],
            num_of_vars: 0,
        }
    }

    fn to_univariate(&self) -> Result<UnivariatePolynomial<F>, String> {
        let mut x_values = vec![];
        for i in 0..self.val.len() {
            x_values.push(F::from(i as u64));
        }
        Ok(UnivariatePolynomial::interpolate(&x_values, &self.val))
    }
}

impl<F: PrimeField> Add for MLE<F> {
    type Output = MLE<F>;

    fn add(self, rhs: Self) -> Self::Output {
        assert!(
            self.val.len() == rhs.val.len(),
            "lhs and rhs must have the same number of evaluations"
        );

        assert!(
            self.val.len().is_power_of_two(),
            "Number of evaluations must be a power of two"
        );

        let mut res = vec![];

        for i in 0..self.val.len() {
            res.push(self.val[i] + rhs.val[i]);
        }

        Self::new(&res)
    }
}

impl<F: PrimeField> Mul for MLE<F> {
    type Output = MLE<F>;

    fn mul(self, rhs: Self) -> Self::Output {
        assert!(
            self.val.len() == rhs.val.len(),
            "lhs and rhs must have the same number of evaluations"
        );

        assert!(
            self.val.len().is_power_of_two(),
            "Number of evaluations must be a power of two"
        );

        let mut res = vec![];

        for i in 0..self.val.len() {
            res.push(self.val[i] * rhs.val[i]);
        }

        Self::new(&res)
    }
}

#[cfg(test)]
mod tests {
    use ark_ff::{Fp64, MontBackend, MontConfig};

    use crate::multilinear_polynomial::{
        coef_form::MultilinearPolynomial, eval_form::MLE, traits::MultilinearPolynomialTrait,
    };

    #[derive(MontConfig)]
    #[modulus = "17"]
    #[generator = "3"]
    pub struct FqConfig;
    pub type Fq = Fp64<MontBackend<FqConfig, 1>>;

    #[test]
    pub fn test_partial_eval_eval_form() {
        let val = vec![
            Fq::from(1),
            Fq::from(2),
            Fq::from(3),
            Fq::from(4),
            Fq::from(5),
            Fq::from(6),
            Fq::from(7),
            Fq::from(8),
        ];

        let poly = MLE::new(&val);

        // number of variables for poly should be 3
        assert!(poly.num_of_vars == 3, "Number of vars for poly should be 3");

        // partially evaluate poly at b = 2
        let reduced_poly = poly.partial_eval(&vec![(2, Fq::from(2))]);

        assert!(reduced_poly.num_of_vars == 2, "Number of vars should be 2");
        assert!(
            dbg!(reduced_poly.val) == vec![Fq::from(5), Fq::from(6), Fq::from(9), Fq::from(10),],
            "Incorrect evaluation"
        );
    }

    #[test]
    fn test_sum_over_boolean_hypercube() {
        let val = vec![
            Fq::from(1),
            Fq::from(2),
            Fq::from(3),
            Fq::from(4),
            Fq::from(5),
            Fq::from(6),
            Fq::from(7),
            Fq::from(8),
        ];

        let poly = MLE::new(&val);

        let res = poly.sum_over_the_boolean_hypercube();

        assert!(
            res == Fq::from(36),
            "Incorrect sum over the boolean hypercube"
        );
    }

    #[test]
    fn test_eval_form_addition() {
        let val1 = vec![
            Fq::from(9),
            Fq::from(12),
            Fq::from(3),
            Fq::from(15),
            Fq::from(24),
            Fq::from(1),
            Fq::from(7),
            Fq::from(9),
        ];

        let val2 = vec![
            Fq::from(1),
            Fq::from(2),
            Fq::from(3),
            Fq::from(4),
            Fq::from(5),
            Fq::from(6),
            Fq::from(7),
            Fq::from(8),
        ];

        let poly1: MLE<ark_ff::Fp<MontBackend<FqConfig, 1>, 1>> = MLE::new(&val1);
        let poly2 = MLE::new(&val2);

        let res_poly = poly1.clone() + poly2.clone();

        let coeff_form = MultilinearPolynomial::interpolate(&res_poly.val);
        assert!(
            poly1.val[3] + poly2.val[3]
                == coeff_form.evaluate(&vec![(0, Fq::from(0)), (1, Fq::from(1)), (2, Fq::from(1))]),
            "Evaluations do not match"
        );
    }

    #[test]
    fn test_eval_form_multiplication() {
        let val1 = vec![
            Fq::from(9),
            Fq::from(12),
            Fq::from(3),
            Fq::from(15),
            Fq::from(24),
            Fq::from(1),
            Fq::from(7),
            Fq::from(9),
        ];

        let val2 = vec![
            Fq::from(1),
            Fq::from(2),
            Fq::from(3),
            Fq::from(4),
            Fq::from(5),
            Fq::from(6),
            Fq::from(7),
            Fq::from(8),
        ];

        let poly1: MLE<ark_ff::Fp<MontBackend<FqConfig, 1>, 1>> = MLE::new(&val1);
        let poly2 = MLE::new(&val2);

        let res_poly = poly1.clone() * poly2.clone();

        let coeff_form = MultilinearPolynomial::interpolate(&res_poly.val);

        assert!(
            poly1.val[3] * poly2.val[3]
                == coeff_form.evaluate(&vec![(0, Fq::from(0)), (1, Fq::from(1)), (2, Fq::from(1))]),
            "Evaluations do not match"
        );
    }

    #[test]
    pub fn test_evaluate_eval_form() {
        let val = vec![
            Fq::from(1),
            Fq::from(2),
            Fq::from(3),
            Fq::from(4),
            Fq::from(5),
            Fq::from(6),
            Fq::from(7),
            Fq::from(8),
        ];

        let poly = MLE::new(&val);

        // number of variables for poly should be 3
        assert!(poly.num_of_vars == 3, "Number of vars for poly should be 3");

        // evaluate poly at a = 3, b = 2 and c = 5
        let res = poly.evaluate(&vec![(1, Fq::from(3)), (3, Fq::from(9)), (2, Fq::from(2))]);

        let coeff_poly = MultilinearPolynomial::interpolate(&val);
        println!("{coeff_poly}");

        assert!(
            dbg!(res)
                == dbg!(coeff_poly.evaluate(&vec![
                    (0, Fq::from(3)),
                    (1, Fq::from(2)),
                    (2, Fq::from(9))
                ])),
            "Incorrect evaluation"
        );
    }

    #[test]
    fn test_eval_form_to_univariate() {
        let evaluations = vec![Fq::from(2), Fq::from(3), Fq::from(8)];
        let eval_poly = MLE::new(&evaluations);
        let eval_poly_univariate = eval_poly.to_univariate().unwrap();

        assert!(
            eval_poly_univariate.evaluate(Fq::from(0)) == Fq::from(2),
            "Incorrect evaluation: Conversion failed"
        );
        assert!(
            eval_poly_univariate.evaluate(Fq::from(1)) == Fq::from(3),
            "Incorrect evaluation: Conversion failed"
        );
        assert!(
            eval_poly_univariate.evaluate(Fq::from(2)) == Fq::from(8),
            "Incorrect evaluation: Conversion failed"
        );
    }
}
