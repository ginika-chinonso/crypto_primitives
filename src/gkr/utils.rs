use ark_ff::PrimeField;

use crate::polynomials::{
    multilinear_poly::{MultilinearPolynomial, MultilinearPolynomialTrait},
    univariate_poly::UnivariatePolynomial,
};

pub fn l_function<F: PrimeField>(b: &[F], c: &[F]) -> Vec<UnivariatePolynomial<F>> {
    b.into_iter()
        .zip(c)
        .map(|(bi, ci)| {
            let m = *ci - *bi;
            // b is the intercept
            // eqn of a straight line = mx + c
            UnivariatePolynomial::<F>::new(vec![*bi, m])
        })
        .collect()
}

pub fn eval_l<F: PrimeField>(l: &[UnivariatePolynomial<F>], r: F) -> Vec<F>{
    l.iter().map(|poly| poly.evaluate(r)).collect()
}

pub fn q_function<F: PrimeField>(
    l: &[UnivariatePolynomial<F>],
    w: &MultilinearPolynomial<F>,
) -> Result<UnivariatePolynomial<F>, &'static str> {
    if w.number_of_vars() != l.len() {
        return Err("w vars and l length should match");
    };

    // Replace each variable for each term in w with its l equivalent polynomial
    // Eg: 2ab where a = l(a) and b = l(b)
    // where l is univariate

    let mut res = UnivariatePolynomial::additive_identity();

    for term in &w.terms {
        let mut coeff = UnivariatePolynomial::new(vec![term.coefficient]);
        for i in 0..term.vars.len() {
            if term.vars[i] {
                coeff = coeff.clone() * l[i].clone();
            }
        }
        res = res + coeff;
    }

    Ok(res)
}

#[cfg(test)]
mod test {

    use crate::{polynomials::{
        multilinear_poly::{
            MultilinearMonomial, MultilinearPolynomial, MultilinearPolynomialTrait,
        },
        univariate_poly::UnivariatePolynomial,
    }, gkr::utils::eval_l};
    use ark_ff::{Fp64, MontBackend, MontConfig};

    use super::{l_function, q_function};

    #[derive(MontConfig)]
    #[modulus = "17"]
    #[generator = "3"]
    pub struct FqConfig;
    pub type Fq = Fp64<MontBackend<FqConfig, 1>>;

    #[test]
    fn test_l_function() {
        let b = vec![Fq::from(2), Fq::from(4), Fq::from(6)];
        let c = vec![Fq::from(3), Fq::from(2), Fq::from(9)];

        let res = l_function(&b, &c);

        assert!(
            res == vec![
                UnivariatePolynomial::new(vec![Fq::from(2), Fq::from(1)]),
                UnivariatePolynomial::new(vec![Fq::from(4), Fq::from(-2)]),
                UnivariatePolynomial::new(vec![Fq::from(6), Fq::from(3)]),
            ]
        )
    }

    #[test]
    fn test_q_function() {
        let l = vec![
            UnivariatePolynomial::new(vec![Fq::from(2), Fq::from(1)]),
            UnivariatePolynomial::new(vec![Fq::from(4), Fq::from(-2)]),
            UnivariatePolynomial::new(vec![Fq::from(6), Fq::from(3)]),
        ];

        let w = MultilinearPolynomial::<Fq>::new(vec![
            MultilinearMonomial::new(Fq::from(2), vec![true, true, false]),
            MultilinearMonomial::new(Fq::from(3), vec![false, true, true]),
        ]);

        let q_func = q_function(&l, &w).unwrap();

        // TODO: Evaluate at a random point
        let q_func_eval = q_func.evaluate(Fq::from(8));

        let l_eval: Vec<Fq> = eval_l(&l, Fq::from(8));

        let w_eval_points = l_eval.into_iter().enumerate().collect();

        let w_eval = w.evaluate(&w_eval_points);

        assert!(q_func_eval == w_eval);
    }
}
