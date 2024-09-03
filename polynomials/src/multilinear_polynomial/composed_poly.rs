use ark_ff::PrimeField;

use crate::{
    multilinear_polynomial::{eval_form::MLE, traits::MultilinearPolynomialTrait},
    univariate_polynomial::UnivariatePolynomial,
};

#[derive(Debug, Clone)]
pub struct ComposedPoly<F: PrimeField + From<i32>> {
    pub polys: Vec<MLE<F>>,
    // operation can be addition or multiplication
    // pub op: String
}

impl<F: PrimeField + From<i32>> ComposedPoly<F> {
    pub fn new(polys: Vec<MLE<F>>) -> Self {
        let len = match polys.get(0) {
            Some(poly) => poly.val.len(),
            None => 0,
        };
        for poly in &polys {
            assert_eq!(poly.val.len(), len);
        }
        Self {
            polys,
            // op
        }
    }

    // Returns the evaluation of the polynomial over the boolean hypercube usually in {0,1}^n
    // This value shouldnt be taken as an outright polynomial because it doesnt capture all the information of the polynomial
    // using it as a polynomial will lead to misleading results
    // This is because the product of two polynomials leads to another polynomial of a  higher degree and would need more than
    // the evaluation over the boolean hypercube of {0,1}^n, where n = number of variables to represent it. A set of {0,1,2,...} may
    // be able to propely represent the higher degree
    pub fn elementwise_mul(&self) -> Vec<F> {
        if self.polys.len() == 0 {
            return vec![];
        } else if self.polys.len() == 1 {
            return self.polys[0].val.to_vec();
        }
        self.polys[0]
            .element_wise_mul(&self.polys[1..].to_vec())
            .val
    }

    pub fn skip_one_and_sum_over_the_boolean_hypercube(&self) -> UnivariatePolynomial<F> {
        let mut res = vec![];

        for i in 0..=self.polys.len() {
            res.push(
                self.partial_eval(&vec![(1, F::from(i as i32))])
                    .elementwise_mul()
                    .iter()
                    .sum(),
            );
        }

        UnivariatePolynomial::interpolate(
            &(0..=self.polys.len())
                .map(|val| F::from(val as i32))
                .collect::<Vec<F>>(),
            &res,
        )
    }
}

impl<F: PrimeField + From<i32>> MultilinearPolynomialTrait<F> for ComposedPoly<F> {
    fn partial_eval(&self, x: &Vec<(usize, F)>) -> Self {
        let partially_evaluated = self.polys.iter().map(|poly| poly.partial_eval(x)).collect();
        Self {
            polys: partially_evaluated,
        }
    }

    fn evaluate(&self, x: &Vec<(usize, F)>) -> F {
        self.polys
            .iter()
            .map(|poly| poly.evaluate(x))
            .fold(F::one(), |init, val| init * val)
    }

    fn number_of_vars(&self) -> usize {
        match self.polys.get(0) {
            Some(poly) => poly.number_of_vars(),
            None => 0,
        }
    }

    fn to_bytes(&self) -> Vec<u8> {
        self.polys.iter().fold(vec![], |mut init, poly| {
            init.append(&mut poly.to_bytes());
            init
        })
    }

    fn relabel(&self) -> Self {
        todo!()
    }

    fn additive_identity() -> Self {
        todo!()
    }

    fn sum_over_the_boolean_hypercube(&self) -> F {
        self.elementwise_mul().iter().sum()
    }

    fn to_univariate(&self) -> Result<UnivariatePolynomial<F>, String> {
        todo!()
        // let mut res = vec![];

        // for i in 0..self.polys.len() {
        //     res.push(self.partial_eval(&vec![(1, F::from(i as i32))]).elementwise_mul().iter().sum());
        // }

        // Ok(UnivariatePolynomial::interpolate(&(0..self.polys.len()).map(|val| F::from(val as i32)).collect::<Vec<F>>(), &res))
    }
}

#[cfg(test)]
pub mod test {
    use ark_bn254::Fq;

    use crate::{
        multilinear_polynomial::{eval_form::MLE, traits::MultilinearPolynomialTrait},
        univariate_polynomial::UnivariatePolynomial,
    };

    use super::ComposedPoly;

    pub fn create_composed_poly() -> ComposedPoly<Fq> {
        // 2ab * 3ab * 5ab = 30a^3b^3
        let poly_1 = MLE::new(&vec![Fq::from(0), Fq::from(0), Fq::from(0), Fq::from(2)]);
        let poly_2 = MLE::new(&vec![Fq::from(0), Fq::from(0), Fq::from(0), Fq::from(3)]);
        let poly_3 = MLE::new(&vec![Fq::from(0), Fq::from(0), Fq::from(0), Fq::from(5)]);

        ComposedPoly::new(vec![poly_1, poly_2, poly_3])
    }

    #[test]
    pub fn test_composed_poly_sum_over_the_boolean_hypercube() {
        let composed_poly = create_composed_poly();

        assert_eq!(
            composed_poly.sum_over_the_boolean_hypercube(),
            Fq::from(30),
            "Invalid sum over the boolean hypercube"
        );
    }

    #[test]
    pub fn test_composed_poly_evaluate() {
        let composed_poly = create_composed_poly();

        let res = composed_poly.evaluate(&vec![(1, Fq::from(3)), (2, Fq::from(5))]);

        assert_eq!(res, Fq::from(101250), "Invalid polynomial evaluation");
    }

    #[test]
    pub fn test_composed_poly_partial_evaluate() {
        let composed_poly = create_composed_poly();

        let res = composed_poly.partial_eval(&vec![(2, Fq::from(5))]);

        assert_eq!(
            res.polys[0].val,
            vec![Fq::from(0), Fq::from(10)],
            "Incorrectpartial evaluation result at index 0"
        );
        assert_eq!(
            res.polys[1].val,
            vec![Fq::from(0), Fq::from(15)],
            "Incorrectpartial evaluation result at index 1"
        );
        assert_eq!(
            res.polys[2].val,
            vec![Fq::from(0), Fq::from(25)],
            "Incorrectpartial evaluation result at index 2"
        );
    }

    #[test]
    pub fn test_composed_poly_skip_one_and_sum_over_the_boolean_hypercube() {
        let composed_poly = create_composed_poly();

        let res = composed_poly.skip_one_and_sum_over_the_boolean_hypercube();

        assert_eq!(
            res,
            UnivariatePolynomial::new(vec![Fq::from(0), Fq::from(0), Fq::from(0), Fq::from(30)]),
            "Incorrect univariate polynomial"
        );
    }
}
