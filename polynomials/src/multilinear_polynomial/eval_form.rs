use crate::{
    multilinear_polynomial::traits::MultilinearPolynomialTrait,
    univariate_polynomial::UnivariatePolynomial, utils::get_sib,
};
use ark_ff::{BigInteger, PrimeField};
use ark_serialize::*;
use std::ops::{Add, Mul};

#[derive(Debug, Clone, CanonicalSerialize, CanonicalDeserialize)]
pub struct MLE<F: PrimeField + From<i32>> {
    // Variables are not zero indexed
    pub num_of_vars: usize,
    // The val vector contains the evaluation of the mle over the boolean hypercube
    pub val: Vec<F>,
}

impl<F: PrimeField + From<i32>> MLE<F> {
    pub fn new(val: &Vec<F>) -> Self {
        assert!(
            val.len().is_power_of_two(),
            "Number of evaluations should be a power of two"
        );
        let num_of_vars = (val.len() as f64).log2().ceil() as usize;
        Self {
            num_of_vars,
            val: val.to_vec(),
        }
    }

    // Indexes are not zero indexed
    pub fn add_variable_at_index(&self, indexes: &mut Vec<usize>) -> MLE<F> {
        if indexes.is_empty() {
            return self.clone();
        }

        indexes.sort();

        let mut new_self = self.val.clone();

        for i in 0..indexes.len() {
            let mut res = vec![F::zero(); new_self.len() * 2];

            let mut shift = 2_usize.pow(((self.num_of_vars + i + 1) - indexes[i]) as u32);
            let mut index = 0;

            for j in 0..new_self.len() {
                if shift == 0 {
                    shift = 2_usize.pow(((self.num_of_vars + i + 1) - indexes[i]) as u32);
                    index += shift;
                }

                res[index] = new_self[j];
                res[get_sib(index, self.num_of_vars + i + 1, indexes[i])] = new_self[j];

                shift -= 1;
                index += 1;
            }

            new_self = res.clone();
        }

        MLE::new(&new_self)
    }

    pub fn skip_one_and_sum_over_the_boolean_hypercube(&self) -> MLE<F> {
        let val1 = self.val[0..self.val.len() / 2]
            .iter()
            .fold(F::zero(), |init, val| init + val);
        let val2 = self.val[(self.val.len() / 2)..]
            .iter()
            .fold(F::zero(), |init, val| init + val);

        Self::new(&vec![val1, val2])
    }

    pub fn skip_one_and_sum_over_the_boolean_hypercube_with_2(&self) -> UnivariatePolynomial<F> {
        let val1 = self.val[0..self.val.len() / 2]
            .iter()
            .fold(F::zero(), |init, val| init + val);
        let val2 = self.val[(self.val.len() / 2)..]
            .iter()
            .fold(F::zero(), |init, val| init + val);
        // (1 - t)y0 + ty1 where t = 2 ; 1 - 2 = -1
        // -1 * val1 = - val1
        let val3 = (F::from(2) * val2) - val1;

        UnivariatePolynomial::interpolate(
            &vec![F::from(0), F::from(1), F::from(2)],
            &vec![val1, val2, val3],
        )
    }

    // This function generates the evaluation form of a polynomial
    // that checks if the values passed in is equal to the value the polynomial was generated with
    // ie, the MLE evaluates to 1 at g and 0 at other places over the boolean hypercube.
    pub fn eq(g: &[F]) -> MLE<F> {
        let num_of_vars = g.len();

        let mut res = vec![F::zero(); 1 << num_of_vars];

        for i in 0..res.len() {
            let binary_string: Vec<F> = format!("{:0width$b}", i, width = num_of_vars)
                .chars()
                .enumerate()
                .map(|(index, bit)| {
                    if bit == '0' {
                        return F::one() - g[index];
                    } else {
                        return g[index];
                    }
                })
                .collect();

            res[i] = binary_string
                .iter()
                .skip(1)
                .fold(binary_string[0], |init, check| init * check);
        }

        MLE::new(&res)
    }

    pub fn skip_one_and_sum_product_over_the_boolean_hypercube_with_2(
        &self,
        rhs: &Self,
    ) -> UnivariatePolynomial<F> {
        let lhs_at_0 = self.val[..self.val.len() / 2].iter();
        let rhs_at_0 = rhs.val[..rhs.val.len() / 2].iter();
        let res_at_0 = lhs_at_0
            .clone()
            .zip(rhs_at_0.clone())
            .fold(F::zero(), |init, (lhs, rhs)| init + (*lhs * *rhs));

        let lhs_at_1 = self.val[self.val.len() / 2..].iter();
        let rhs_at_1 = rhs.val[rhs.val.len() / 2..].iter();
        let res_at_1 = lhs_at_1
            .clone()
            .zip(rhs_at_1.clone())
            .fold(F::zero(), |init, (lhs, rhs)| init + (*lhs * *rhs));

        let lhs_at_2 = lhs_at_0
            .zip(lhs_at_1)
            .map(|(lhs, rhs)| (F::one().neg() * lhs) + (F::from(2) * rhs));
        let rhs_at_2 = rhs_at_0
            .zip(rhs_at_1)
            .map(|(lhs, rhs)| (F::one().neg() * lhs) + (F::from(2) * rhs));
        let res_at_2 = lhs_at_2
            .zip(rhs_at_2)
            .fold(F::zero(), |init, (lhs, rhs)| init + (lhs * rhs));

        UnivariatePolynomial::interpolate(
            &vec![F::from(0), F::from(1), F::from(2)],
            &vec![res_at_0, res_at_1, res_at_2],
        )
    }

    pub fn sum_product_over_the_boolean_hypercube(&self, rhs: &MLE<F>) -> F {
        self.val
            .iter()
            .zip(&rhs.val)
            .fold(F::zero(), |init, (lhs, rhs)| init + (*lhs * rhs))
    }

    pub fn element_wise_mul(&self, rhs: &Self) -> Self {
        assert_eq!(
            self.val.len(),
            rhs.val.len(),
            "LHS and RHS should have same number of evaluations"
        );

        let res_arr: Vec<F> = self
            .val
            .iter()
            .zip(&rhs.val)
            .map(|(lhs, rhs)| *lhs * rhs)
            .collect();

        MLE::new(&res_arr)
    }
}

impl<F: PrimeField + From<i32>> MultilinearPolynomialTrait<F> for MLE<F> {
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
            let mut shift = 2_usize.pow((num_of_vars - var) as u32);

            for _ in 0..(new_poly.val.len() / 2) {
                if shift == 0 {
                    shift = 2_usize.pow((num_of_vars - var) as u32);
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
        assert!(
            points.len() == self.num_of_vars,
            "Provide evaluation point for all variables"
        );
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
        self.val.iter().fold(
            self.num_of_vars.to_be_bytes().to_vec(),
            |mut init, value| {
                init.extend(value.into_bigint().to_bytes_be());
                init
            },
        )
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

impl<F: PrimeField + From<i32>> Add for MLE<F> {
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

// Multiplication takes variabes as different variables
// Eg: ab * ab = abcd
impl<F: PrimeField + From<i32>> Mul for MLE<F> {
    type Output = MLE<F>;

    fn mul(self, rhs: Self) -> Self::Output {
        assert!(
            self.val.len().is_power_of_two() && rhs.val.len().is_power_of_two(),
            "Length of evaluations should be a power of two"
        );

        let res_num_of_vars = self.number_of_vars() + rhs.number_of_vars();

        let mut ind_self = (1..=(res_num_of_vars - self.number_of_vars())).collect::<Vec<_>>();
        let mut ind_rhs = (rhs.number_of_vars() + 1..=res_num_of_vars).collect::<Vec<_>>();

        let new_self = self.add_variable_at_index(&mut ind_self);
        let new_rhs = rhs.add_variable_at_index(&mut ind_rhs);

        let mut res = vec![F::zero(); new_self.val.len()];

        for i in 0..new_self.val.len() {
            res[i] = new_self.val[i] * new_rhs.val[i];
        }

        Self::new(&res)
    }
}

//////////////////////////////////
/// EVAL FORM SPARSE POLY
/// //////////////////////////////
pub mod sparse_mle {
    use ark_ff::PrimeField;

    use crate::{
        multilinear_polynomial::traits::MultilinearPolynomialTrait,
        univariate_polynomial::UnivariatePolynomial,
    };

    #[derive(Debug)]
    pub struct SparseMle<F> {
        pub num_of_variables: usize,
        pub values: Vec<(usize, F)>,
    }

    impl<F: PrimeField> SparseMle<F> {
        pub fn new(num_of_variables: usize, values: Vec<(usize, F)>) -> Self {
            Self {
                num_of_variables,
                values,
            }
        }
    }

    impl<F: PrimeField> MultilinearPolynomialTrait<F> for SparseMle<F> {
        fn partial_eval(&self, x: &Vec<(usize, F)>) -> Self {
            todo!()
        }

        fn evaluate(&self, x: &Vec<(usize, F)>) -> F {
            todo!()
        }

        fn number_of_vars(&self) -> usize {
            self.num_of_variables
        }

        fn to_bytes(&self) -> Vec<u8> {
            todo!()
        }

        fn relabel(&self) -> Self {
            todo!()
        }

        fn additive_identity() -> Self {
            todo!()
        }

        fn sum_over_the_boolean_hypercube(&self) -> F {
            self.values
                .iter()
                .fold(F::zero(), |init, (_, eval)| init + eval)
        }

        fn to_univariate(&self) -> Result<UnivariatePolynomial<F>, String> {
            todo!()
        }
    }

    #[cfg(test)]
    pub mod tests {
        use ark_bn254::Fq;

        use super::{MultilinearPolynomialTrait, SparseMle};

        #[test]
        pub fn test_sparse_poly_sum_over_the_boolean_hypercube() {
            let values = vec![(3, Fq::from(1)), (7, Fq::from(1))];

            let sparse_poly = SparseMle::new(3, values);

            assert_eq!(
                sparse_poly.sum_over_the_boolean_hypercube(),
                Fq::from(2),
                "Incorrect sum over the boolean hypercube"
            );
        }
    }
}

#[cfg(test)]
mod tests {
    use ark_bn254::Fr;

    use crate::multilinear_polynomial::{
        coef_form::MultilinearPolynomial, eval_form::MLE, traits::MultilinearPolynomialTrait,
    };

    pub type Fq = Fr;

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
            reduced_poly.val == vec![Fq::from(5), Fq::from(6), Fq::from(9), Fq::from(10),],
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

        let poly1: MLE<Fq> = MLE::new(&val1);
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
    pub fn test_evaluate_eval_form() {
        // Polynomial in consideration = 2ab + 3bc
        let val = vec![
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(3),
            Fq::from(0),
            Fq::from(0),
            Fq::from(2),
            Fq::from(5),
        ];

        let poly = MLE::new(&val);

        // number of variables for poly should be 3
        assert!(poly.num_of_vars == 3, "Number of vars for poly should be 3");

        // evaluate poly at a = 3, b = 2 and c = 5
        let res = poly.evaluate(&vec![(1, Fq::from(3)), (3, Fq::from(5)), (2, Fq::from(2))]);

        assert!(res == Fq::from(42), "Incorrect evaluation");
    }

    #[test]
    fn test_eval_form_to_univariate() {
        let evaluations = vec![Fq::from(2), Fq::from(3), Fq::from(8), Fq::from(12)];
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
        assert!(
            eval_poly_univariate.evaluate(Fq::from(3)) == Fq::from(12),
            "Incorrect evaluation: Conversion failed"
        );
    }

    #[test]
    fn test_eval_form_multiplication() {
        // let poly1 = 2ab
        let val1 = vec![Fq::from(0), Fq::from(0), Fq::from(0), Fq::from(2)];

        // let poly2 = 3cd
        let val2 = vec![Fq::from(0), Fq::from(0), Fq::from(0), Fq::from(3)];

        let poly1: MLE<Fq> = MLE::new(&val1);
        let poly2: MLE<Fq> = MLE::new(&val2);

        // the resulting poly should be 6abcd
        let res_poly = poly1.clone() * poly2.clone();

        assert!(
            poly1.evaluate(&vec![
                (1, Fq::from(38)),
                (2, Fq::from(64)),
                // (3, Fq::from(90)),
                // (4, Fq::from(30))
            ]) * poly2.evaluate(&vec![
                // (1, Fq::from(38)),
                // (2, Fq::from(64)),
                (1, Fq::from(90)),
                (2, Fq::from(30))
            ]) == res_poly.evaluate(&vec![
                (1, Fq::from(38)),
                (2, Fq::from(64)),
                (3, Fq::from(90)),
                (4, Fq::from(30))
            ]),
            "Evaluations do not match"
        );
    }

    #[test]
    pub fn test_add_variable_at_index_1() {
        // Polynomial in consideration: 2ab
        let val = vec![Fq::from(0), Fq::from(0), Fq::from(0), Fq::from(2)];

        let poly = MLE::new(&val);

        // add a new variable to the front to get 2abc where a is the new variable

        // Note that indexes are not zero indexed
        let new_poly = poly.add_variable_at_index(&mut vec![1]);

        assert!(
            new_poly.val
                == vec![
                    Fq::from(0),
                    Fq::from(0),
                    Fq::from(0),
                    Fq::from(2),
                    Fq::from(0),
                    Fq::from(0),
                    Fq::from(0),
                    Fq::from(2),
                ],
            "Failed to add variable"
        );
    }

    #[test]
    pub fn test_add_variable_at_index_2() {
        // Polynomial in consideration: 2ab
        let val = vec![Fq::from(0), Fq::from(0), Fq::from(0), Fq::from(2)];

        let poly = MLE::new(&val);

        // add a new variable at the middle to get 2abc where b is the new variable

        // Note that indexes are not zero indexed
        let new_poly = poly.add_variable_at_index(&mut vec![2]);

        assert!(
            new_poly.val
                == vec![
                    Fq::from(0),
                    Fq::from(0),
                    Fq::from(0),
                    Fq::from(0),
                    Fq::from(0),
                    Fq::from(2),
                    Fq::from(0),
                    Fq::from(2),
                ],
            "Failed to add variable"
        );
    }

    #[test]
    pub fn test_add_variable_at_index_3() {
        // Polynomial in consideration: 2ab
        let val = vec![Fq::from(0), Fq::from(0), Fq::from(0), Fq::from(2)];

        let poly = MLE::new(&val);

        // add a new variable at the end to get 2abc where c is the new variable

        // Note that indexes are not zero indexed
        let new_poly = poly.add_variable_at_index(&mut vec![3]);

        assert!(
            new_poly.val
                == vec![
                    Fq::from(0),
                    Fq::from(0),
                    Fq::from(0),
                    Fq::from(0),
                    Fq::from(0),
                    Fq::from(0),
                    Fq::from(2),
                    Fq::from(2),
                ],
            "Failed to add variable"
        );
    }

    #[test]
    pub fn test_add_variable_at_index_1_and_2() {
        // Polynomial in consideration: 2ab
        let val = vec![Fq::from(0), Fq::from(0), Fq::from(0), Fq::from(2)];

        let poly = MLE::new(&val);

        // add a new variable at index 1 and 2 to get 2abcd where a and b are the new variable

        // Note that indexes are not zero indexed
        let new_poly = poly.add_variable_at_index(&mut vec![2, 1]);

        assert!(
            new_poly.val
                == vec![
                    Fq::from(0),
                    Fq::from(0),
                    Fq::from(0),
                    Fq::from(2),
                    Fq::from(0),
                    Fq::from(0),
                    Fq::from(0),
                    Fq::from(2),
                    Fq::from(0),
                    Fq::from(0),
                    Fq::from(0),
                    Fq::from(2),
                    Fq::from(0),
                    Fq::from(0),
                    Fq::from(0),
                    Fq::from(2),
                ],
            "Failed to add variable"
        );
    }

    #[test]
    pub fn test_add_variable_at_index_2_and_3() {
        // Polynomial in consideration: 2ab
        let val = vec![Fq::from(0), Fq::from(0), Fq::from(0), Fq::from(2)];

        let poly = MLE::new(&val);

        // add a new variable at index 2 and 3 to get 2abcd where b and c are the new variables

        // Note that indexes are not zero indexed
        let new_poly = poly.add_variable_at_index(&mut vec![3, 2]);

        assert!(
            new_poly.val
                == vec![
                    Fq::from(0),
                    Fq::from(0),
                    Fq::from(0),
                    Fq::from(0),
                    Fq::from(0),
                    Fq::from(0),
                    Fq::from(0),
                    Fq::from(0),
                    Fq::from(0),
                    Fq::from(2),
                    Fq::from(0),
                    Fq::from(2),
                    Fq::from(0),
                    Fq::from(2),
                    Fq::from(0),
                    Fq::from(2),
                ],
            "Failed to add variable"
        );
    }

    #[test]
    pub fn test_add_variable_at_index_3_and_4() {
        // Polynomial in consideration: 2ab
        let val = vec![Fq::from(0), Fq::from(0), Fq::from(0), Fq::from(2)];

        let poly = MLE::new(&val);

        // add a new variable at index 3 and 4 to get 2abcd where c and d are the new variables

        // Note that indexes are not zero indexed
        let new_poly = poly.add_variable_at_index(&mut vec![4, 3]);

        assert!(
            new_poly.val
                == vec![
                    Fq::from(0),
                    Fq::from(0),
                    Fq::from(0),
                    Fq::from(0),
                    Fq::from(0),
                    Fq::from(0),
                    Fq::from(0),
                    Fq::from(0),
                    Fq::from(0),
                    Fq::from(0),
                    Fq::from(0),
                    Fq::from(0),
                    Fq::from(2),
                    Fq::from(2),
                    Fq::from(2),
                    Fq::from(2),
                ],
            "Failed to add variable"
        );
    }

    #[test]
    pub fn test_add_variable_at_index_1_and_4() {
        // Polynomial in consideration: 2ab
        let val = vec![Fq::from(0), Fq::from(0), Fq::from(0), Fq::from(2)];

        let poly = MLE::new(&val);

        // add a new variable at index 1 and 4 to get 2abcd where a and d are the new variables

        // Note that indexes are not zero indexed
        let new_poly = poly.add_variable_at_index(&mut vec![4, 1]);

        assert!(
            new_poly.val
                == vec![
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
                ],
            "Failed to add variable"
        );
    }

    #[test]
    pub fn test_skip_one_and_sum_over_the_boolean_hypercube() {
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

        let poly = MLE::new(&val);

        let claimed_sum = poly.sum_over_the_boolean_hypercube();

        let univariate_poly = poly.skip_one_and_sum_over_the_boolean_hypercube();

        assert!(
            univariate_poly.val == vec![Fq::from(4), Fq::from(4)],
            "Wrong univariate polynomial"
        );

        assert!(
            claimed_sum
                == univariate_poly.evaluate(&vec![(1, Fq::from(0))])
                    + univariate_poly.evaluate(&vec![(1, (Fq::from(1)))]),
            "Invalid univariate poly"
        );
    }

    #[test]
    pub fn test_element_wise_mul() {
        let vec_1 = vec![Fq::from(0), Fq::from(0), Fq::from(0), Fq::from(2)];

        let vec_2 = vec![Fq::from(0), Fq::from(0), Fq::from(0), Fq::from(3)];

        let poly_1 = MLE::new(&vec_1);

        let poly_2 = MLE::new(&vec_2);

        let res_poly = poly_1.element_wise_mul(&poly_2);

        assert_eq!(
            res_poly.val,
            vec![Fq::from(0), Fq::from(0), Fq::from(0), Fq::from(6),],
            "Wrong elementwise mul result"
        );
    }

    #[test]
    pub fn test_skip_one_and_sum_product_over_the_boolean_hypercube_with_2() {
        let val_1 = vec![Fq::from(0), Fq::from(0), Fq::from(0), Fq::from(2)];

        let val_2 = vec![Fq::from(0), Fq::from(0), Fq::from(0), Fq::from(3)];

        let poly_1 = MLE::new(&val_1);
        let poly_2 = MLE::new(&val_2);

        let res_poly = poly_1.skip_one_and_sum_product_over_the_boolean_hypercube_with_2(&poly_2);

        assert_eq!(
            res_poly.evaluate(Fq::from(0)),
            Fq::from(0),
            "Invalid evaluation at 0 of sum of product over the boolean hypercube"
        );
        assert_eq!(
            res_poly.evaluate(Fq::from(1)),
            Fq::from(6),
            "Invalid evaluation at 1 of sum of product over the boolean hypercube"
        );
        assert_eq!(
            res_poly.evaluate(Fq::from(2)),
            Fq::from(24),
            "Invalid evaluation at 2 sum of product over the boolean hypercube"
        );
    }

    #[test]
    pub fn test_eq_function() {
        let g = [Fq::from(1), Fq::from(0), Fq::from(1)];

        let eq_poly = MLE::eq(&g);

        dbg!(&eq_poly.num_of_vars);

        let sum = eq_poly.sum_over_the_boolean_hypercube();

        assert_eq!(sum, Fq::from(1), "Incorrect sum over the boolean hypercube");

        assert_eq!(
            Fq::from(1),
            eq_poly.evaluate(
                &g.into_iter()
                    .enumerate()
                    .map(|(var, val)| (var + 1, val))
                    .collect::<Vec<(usize, Fq)>>()
            ),
            "Eq poly evaluated at g shoild give 1"
        );
    }

    #[test]
    pub fn test_skip_one_and_sum_over_the_boolean_hypercube_with_2() {
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

        let poly = MLE::new(&val);

        let univariate_poly = poly.skip_one_and_sum_over_the_boolean_hypercube_with_2();

        assert_eq!(
            univariate_poly.evaluate(Fq::from(0)),
            Fq::from(4),
            "Wrong univariate polynomial evaluation"
        );
        assert_eq!(
            univariate_poly.evaluate(Fq::from(0)),
            Fq::from(4),
            "Wrong univariate polynomial evaluation"
        );
        assert_eq!(
            univariate_poly.evaluate(Fq::from(0)),
            Fq::from(4),
            "Wrong univariate polynomial evaluation"
        );
    }

    #[test]
    pub fn test_sum_product_over_the_boolean_hypercube_with_2() {
        let poly1 = MLE::new(&vec![Fq::from(0), Fq::from(0), Fq::from(0), Fq::from(2)]);
        let poly2 = MLE::new(&vec![Fq::from(0), Fq::from(0), Fq::from(0), Fq::from(3)]);

        let claimed_sum = poly1.sum_product_over_the_boolean_hypercube(&poly2);

        assert_eq!(
            claimed_sum,
            Fq::from(6),
            "Incorrect sum of product over the boolean hypercube"
        );
    }
}

// type MLE = Vec<u64>;

// enum Polynomial {
//     PMLE(MLE),
//     Product([MLE; 2])
// }

// constraint
// product polynomials should be uniform (i.e the same number of variables)
// and each variable should match

// composed evaluate
// self.mles.fold(|acc, mle| acc *= mle.eval(point), 1)

// partial -> Polynomial
// self.mles.map(|mle| mle.partial_eval(point)).into()

// what does sumcheck need?
// n_vars :ok
// for every var, generate round poly (deg = 2)
//  3 points => 0_half, 1_half and 2_half
//. 0_half => partial eval at a = 0
//. 1_half => partial eval at a = 1
//. 2_half => partial eval at a = 2

// example
// 00 00
// 01 01
// 10 10
// 11 11

// 0_half
// [00, 01] hadamard_prod [00, 01]

// 1_half
// 10, 11

// 2_half
// 20, 21

// round_poly = [0_half.sum(), 1_half.sum(), 2_half.sum()]

// let mut round_poly = vec![];
// for i in 0..(deg+1) {
//  round_poly.push(poly.get_half(i).sum())
//}

// fn main() {
//     dbg!(0 << 1);
// }
