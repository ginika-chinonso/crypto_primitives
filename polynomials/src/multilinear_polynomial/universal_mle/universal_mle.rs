use ark_ff::PrimeField;

use crate::univariate_polynomial::UnivariatePolynomial;

use super::{
    super::{eval_form::MLE, traits::MultilinearPolynomialTrait},
    ops::Op,
};

#[derive(Debug, Clone)]
pub enum UniversalMLE<F: PrimeField + From<i32>> {
    MLES(Vec<UniversalMLE<F>>, Op<F>),
    Value(MLE<F>),
}

impl<F: PrimeField + From<i32>> UniversalMLE<F> {
    pub fn new(mles: Vec<UniversalMLE<F>>, op: Op<F>) -> UniversalMLE<F> {
        // TODO: verify that all mles have the same number of variables
        mles.iter().skip(1).fold(
            mles.get(0)
                .unwrap_or(&UniversalMLE::additive_identity())
                .number_of_vars(),
            |init, mle| {
                assert_eq!(
                    init,
                    mle.number_of_vars(),
                    "Number of variables do not match"
                );
                init
            },
        );
        UniversalMLE::MLES(mles, op)
    }

    pub fn to_eval(&self) -> Vec<F> {
        match self {
            UniversalMLE::MLES(mles, op) => op.reduce(
                &mles
                    .iter()
                    .map(|mle| mle.to_eval())
                    .collect::<Vec<Vec<F>>>(),
            ),
            UniversalMLE::Value(mle) => mle.val.to_vec(),
        }
    }

    pub fn skip_one_and_sum_over_the_boolean_hypercube(&self) -> UnivariatePolynomial<F> {
        let res: Vec<Vec<F>> = match self {
            UniversalMLE::MLES(mles, op) => (0..=(op.degree_checker)(mles))
                .map(|num| match F::from_str(&num.to_string()) {
                    Ok(ind) => self.partial_eval(&vec![(1, ind)]).to_eval(),
                    Err(_) => panic!("Error converting value to string"),
                })
                .collect(),
            UniversalMLE::Value(mle) => (0..2)
                .map(|ind| match F::from_str(&ind.to_string()) {
                    Ok(val) => mle.partial_eval(&vec![(1, val)]).val,
                    Err(_) => panic!("Error converting number to string"),
                })
                .collect(),
        };
        let y_values = res.iter().map(|arr| arr.iter().sum()).collect();
        let x_values = (0..=self.degree())
            .map(|ind| match F::from_str(&ind.to_string()) {
                Ok(val) => val,
                Err(_) => panic!("Error converting number to string"),
            })
            .collect();

        UnivariatePolynomial::interpolate(&x_values, &y_values)
    }

    pub fn degree(&self) -> usize {
        match self {
            UniversalMLE::MLES(mles, op) => (op.degree_checker)(mles),
            UniversalMLE::Value(_) => 1,
        }
    }
}

impl<F: PrimeField + From<i32>> MultilinearPolynomialTrait<F> for UniversalMLE<F> {
    fn partial_eval(&self, x: &Vec<(usize, F)>) -> UniversalMLE<F> {
        match self {
            UniversalMLE::MLES(mles, op) => {
                let partially_evaluated_mles = mles
                    .iter()
                    .map(|mle| mle.partial_eval(x))
                    .collect::<Vec<UniversalMLE<F>>>();
                UniversalMLE::MLES(partially_evaluated_mles, op.clone())
            }
            UniversalMLE::Value(mle) => {
                let res = mle.partial_eval(x);
                UniversalMLE::Value(res)
            }
        }
    }

    fn evaluate(&self, x: &Vec<(usize, F)>) -> F {
        match self {
            UniversalMLE::MLES(mles, op) => {
                let evaluated_results = mles.iter().map(|mle| mle.evaluate(x)).collect::<Vec<F>>();
                let first = evaluated_results[0];
                evaluated_results
                    .iter()
                    .skip(1)
                    .fold(first, |init, evaluation| (op.op)(init, *evaluation))
            }
            UniversalMLE::Value(mle) => mle.evaluate(x),
        }
    }

    fn number_of_vars(&self) -> usize {
        match self {
            UniversalMLE::MLES(mles, _) => match mles.get(0) {
                Some(mle) => mle.number_of_vars(),
                None => 0,
            },
            UniversalMLE::Value(mle) => mle.number_of_vars(),
        }
    }

    fn to_bytes(&self) -> Vec<u8> {
        match self {
            UniversalMLE::MLES(mles, op) => mles.iter().fold(op.to_bytes(), |mut init, mle| {
                init.extend(mle.to_bytes());
                init
            }),
            UniversalMLE::Value(mle) => mle.to_bytes(),
        }
    }

    fn relabel(&self) -> Self {
        todo!()
    }

    fn additive_identity() -> Self {
        UniversalMLE::Value(MLE::new(&vec![F::zero()]))
    }

    fn sum_over_the_boolean_hypercube(&self) -> F {
        match self {
            UniversalMLE::MLES(mles, op) => op
                .reduce(
                    &mles
                        .iter()
                        .map(|mle| mle.to_eval().to_vec())
                        .collect::<Vec<Vec<F>>>(),
                )
                .iter()
                .sum(),
            UniversalMLE::Value(mle) => mle.sum_over_the_boolean_hypercube(),
        }
    }

    fn to_univariate(&self) -> Result<UnivariatePolynomial<F>, String> {
        todo!()
    }
}

#[cfg(test)]
pub mod test {

    use crate::{
        multilinear_polynomial::{
            eval_form::MLE, traits::MultilinearPolynomialTrait, universal_mle::ops::Ops,
        },
        univariate_polynomial::UnivariatePolynomial,
    };
    use ark_bn254::Fq;

    use super::UniversalMLE;

    pub fn create_add_universal_mle() -> UniversalMLE<Fq> {
        // 2ab + 4
        let poly1 = UniversalMLE::Value(MLE::new(&vec![
            Fq::from(4),
            Fq::from(4),
            Fq::from(4),
            Fq::from(6),
        ]));
        // 3ab
        let poly2 = UniversalMLE::Value(MLE::new(&vec![
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(3),
        ]));
        // 5ab - 1
        let poly3 = UniversalMLE::Value(MLE::new(&vec![
            Fq::from(-1),
            Fq::from(-1),
            Fq::from(-1),
            Fq::from(4),
        ]));

        UniversalMLE::new(vec![poly1, poly2, poly3], Ops::new().add)
    }

    pub fn create_mul_universal_mle() -> UniversalMLE<Fq> {
        // 2ab + 4
        let poly1 = UniversalMLE::Value(MLE::new(&vec![
            Fq::from(4),
            Fq::from(4),
            Fq::from(4),
            Fq::from(6),
        ]));
        // 3ab
        let poly2 = UniversalMLE::Value(MLE::new(&vec![
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(3),
        ]));
        // 5ab - 1
        let poly3 = UniversalMLE::Value(MLE::new(&vec![
            Fq::from(-1),
            Fq::from(-1),
            Fq::from(-1),
            Fq::from(4),
        ]));

        UniversalMLE::new(vec![poly1, poly2, poly3], Ops::new().mul)
    }

    pub fn create_higher_degree_poly() -> UniversalMLE<Fq> {
        let poly1 = UniversalMLE::Value(MLE::new(&vec![
            Fq::from(2),
            Fq::from(2),
            Fq::from(2),
            Fq::from(2),
        ]));

        let poly2 = UniversalMLE::Value(MLE::new(&vec![
            Fq::from(0),
            Fq::from(0),
            Fq::from(3),
            Fq::from(3),
        ]));

        let poly3 = UniversalMLE::Value(MLE::new(&vec![
            Fq::from(0),
            Fq::from(2),
            Fq::from(3),
            Fq::from(5),
        ]));

        let poly4 = UniversalMLE::Value(MLE::new(&vec![
            Fq::from(0),
            Fq::from(0),
            Fq::from(1),
            Fq::from(1),
        ]));

        let poly1_mul_poly2 = UniversalMLE::new(vec![poly1, poly2], Ops::new().mul);
        let poly3_mul_poly4 = UniversalMLE::new(vec![poly3, poly4], Ops::new().mul);

        UniversalMLE::new(
            vec![
                poly1_mul_poly2.clone(),
                poly1_mul_poly2,
                poly3_mul_poly4.clone(),
                poly3_mul_poly4,
            ],
            Ops::new().mul,
        )
    }

    pub fn create_mul_and_add_polys() -> UniversalMLE<Fq> {
        let poly1 = UniversalMLE::Value(MLE::new(&vec![
            Fq::from(2),
            Fq::from(2),
            Fq::from(2),
            Fq::from(2),
            Fq::from(2),
            Fq::from(2),
            Fq::from(2),
            Fq::from(5),
        ]));

        // 2a + 3b
        let poly2 = UniversalMLE::Value(MLE::new(&vec![
            Fq::from(0),
            Fq::from(0),
            Fq::from(3),
            Fq::from(3),
            Fq::from(2),
            Fq::from(2),
            Fq::from(5),
            Fq::from(5),
        ]));

        // 3b + 2c
        let poly3 = UniversalMLE::Value(MLE::new(&vec![
            Fq::from(0),
            Fq::from(2),
            Fq::from(3),
            Fq::from(5),
            Fq::from(0),
            Fq::from(2),
            Fq::from(3),
            Fq::from(5),
        ]));

        // 5ac + b
        let poly4 = UniversalMLE::Value(MLE::new(&vec![
            Fq::from(0),
            Fq::from(0),
            Fq::from(1),
            Fq::from(1),
            Fq::from(0),
            Fq::from(5),
            Fq::from(1),
            Fq::from(6),
        ]));

        let poly1_mul_poly3 = UniversalMLE::new(vec![poly1, poly3], Ops::new().mul);
        let poly2_mul_poly4 = UniversalMLE::new(vec![poly2, poly4], Ops::new().mul);

        UniversalMLE::new(vec![poly1_mul_poly3, poly2_mul_poly4], Ops::new().add)
    }

    #[test]
    pub fn test_add_universalmle_partial_eval() {
        let universal_poly = create_add_universal_mle();
        // partially evaluate at a = 3
        let res = universal_poly.partial_eval(&vec![(1, Fq::from(3))]);
        match res {
            UniversalMLE::MLES(mles, _) => {
                match &mles[0] {
                    UniversalMLE::Value(mle) => {
                        assert_eq!(
                            mle.val,
                            vec![Fq::from(4), Fq::from(10)],
                            "Wrong partial evaluation of 2ab + 4"
                        );
                    }
                    UniversalMLE::MLES(_, _) => {}
                }
                match &mles[1] {
                    UniversalMLE::Value(mle) => {
                        assert_eq!(
                            mle.val,
                            vec![Fq::from(0), Fq::from(9)],
                            "Wrong partial evaluation of 3ab"
                        );
                    }
                    UniversalMLE::MLES(_, _) => {}
                }
                match &mles[2] {
                    UniversalMLE::Value(mle) => {
                        assert_eq!(
                            mle.val,
                            vec![Fq::from(-1), Fq::from(14)],
                            "Wrong partial evaluation of 5ab - 1"
                        );
                    }
                    UniversalMLE::MLES(_, _) => {}
                }
            }
            UniversalMLE::Value(_) => {}
        }
    }

    #[test]
    pub fn test_mul_universalmle_partial_eval() {
        let poly = create_mul_universal_mle();
        // partially evaluate at b = 4
        let res: UniversalMLE<Fq> = poly.partial_eval(&vec![(2, Fq::from(4))]);
        match res {
            UniversalMLE::MLES(mle, _) => {
                match &mle[0] {
                    UniversalMLE::Value(mle) => {
                        assert_eq!(
                            mle.val,
                            vec![Fq::from(4), Fq::from(12)],
                            "Invalid poly partial eval"
                        );
                    }
                    UniversalMLE::MLES(_, _) => {}
                }

                match &mle[1] {
                    UniversalMLE::Value(mle) => {
                        assert_eq!(
                            mle.val,
                            vec![Fq::from(0), Fq::from(12)],
                            "Invalid poly partial evaluation"
                        );
                    }
                    UniversalMLE::MLES(_, _) => {}
                }

                match &mle[2] {
                    UniversalMLE::Value(mle) => {
                        assert_eq!(
                            mle.val,
                            vec![Fq::from(-1), Fq::from(19)],
                            "Invalid poly partial evaluation"
                        );
                    }
                    UniversalMLE::MLES(_, _) => {}
                }
            }
            UniversalMLE::Value(_) => {}
        }
    }

    #[test]
    pub fn test_add_universalmle_eval() {
        let poly = create_add_universal_mle();
        // evaluate at a = 3 and b = 4
        // res = 10ab + 3
        let res = poly.evaluate(&vec![(1, Fq::from(3)), (2, Fq::from(4))]);
        assert_eq!(res, Fq::from(123), "Invalid evaluation");
    }

    #[test]
    pub fn test_mul_universalmle_eval() {
        let poly = create_mul_universal_mle();
        // evaluate at a = 3 and b = 4
        // res = 30a^3b^3 + 60a^2b^2 - 6a^2b^2 - 12ab
        let res = poly.evaluate(&vec![(1, Fq::from(3)), (2, Fq::from(4))]);
        assert_eq!(res, Fq::from(59472), "Invalid evaluation");
    }

    #[test]
    pub fn test_add_universalmle_sum_over_the_boolean_hypercube() {
        let poly = create_add_universal_mle();

        // poly when combined should give 10ab + 3
        let res = poly.sum_over_the_boolean_hypercube();

        assert_eq!(res, Fq::from(22), "Invalid sum over the boolean hypercube");
    }

    #[test]
    pub fn test_mul_universalmle_sum_over_the_boolean_hypercube() {
        let poly = create_mul_universal_mle();

        // res = 30a^3b^3 + 60a^2b^2 - 6a^2b^2 - 12ab
        let res = poly.sum_over_the_boolean_hypercube();

        assert_eq!(res, Fq::from(72), "Invalid sum over the boolean hypercube");
    }

    #[test]
    pub fn test_op_reduce() {
        let poly = create_mul_and_add_polys();

        let reduced_form = poly.to_eval();

        assert_eq!(
            reduced_form,
            vec![
                Fq::from(0),
                Fq::from(4),
                Fq::from(9),
                Fq::from(13),
                Fq::from(0),
                Fq::from(14),
                Fq::from(11),
                Fq::from(55),
            ],
            "Incorrect poly reduction"
        );
    }

    #[test]
    pub fn test_op_degree_checker() {
        let poly = create_mul_and_add_polys();

        assert_eq!(poly.degree(), 2, "Incorrect polynomial degree");
    }

    #[test]
    pub fn test_gkr_like_poly() {
        let poly = create_mul_and_add_polys();

        // res = 9ab^2c + 6abc^2 + 6b + 4c + 10a^2c + 2ab + 15abc + 3b^2
        let res = poly.sum_over_the_boolean_hypercube();
        assert_eq!(
            res,
            Fq::from(106),
            "Incorrect sum over the boolean hypercube"
        );
    }

    #[test]
    pub fn test_universalmle_degree() {
        let poly = create_mul_and_add_polys();

        assert_eq!(poly.degree(), 2, "Incorrect polynomial degree");
    }

    #[test]
    pub fn test_universalmle_higher_degree() {
        assert_eq!(
            create_mul_universal_mle().degree(),
            3,
            "Incorrect polynomial degree"
        );
        assert_eq!(
            create_higher_degree_poly().degree(),
            8,
            "Incorrect polynomial degree"
        );
    }

    #[test]
    pub fn test_universalmle_skip_one_and_sum_over_the_boolean_hypercube() {
        let poly = create_mul_and_add_polys();

        let res = poly.skip_one_and_sum_over_the_boolean_hypercube();

        assert_eq!(
            res,
            UnivariatePolynomial::interpolate(
                &vec![Fq::from(0), Fq::from(1), Fq::from(2)],
                &vec![Fq::from(26), Fq::from(80), Fq::from(174)]
            ),
            "Invalid sum over the boolean hypercube"
        );
    }
}
