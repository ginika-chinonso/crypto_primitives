use std::cmp::max;

use ark_ff::PrimeField;

use super::universal_mle::UniversalMLE;

#[derive(Debug, PartialEq, Clone)]
pub struct Op<F: PrimeField + From<i32>> {
    pub name: String,
    pub op: fn(F, F) -> F,
    pub degree_checker: fn(&Vec<UniversalMLE<F>>) -> usize,
}

impl<F: PrimeField + From<i32>> Op<F> {
    pub fn new(
        name: String,
        op: fn(F, F) -> F,
        degree_checker: fn(&Vec<UniversalMLE<F>>) -> usize,
    ) -> Self {
        Self {
            name,
            op,
            degree_checker,
        }
    }

    pub fn reduce(&self, mles_eval: &Vec<Vec<F>>) -> Vec<F> {
        mles_eval
            .iter()
            .skip(1)
            .fold(mles_eval[0].to_vec(), |mut init, mle_eval| {
                init = init
                    .iter()
                    .zip(mle_eval)
                    .map(|(lhs, rhs)| (self.op)(*lhs, *rhs))
                    .collect();
                init
            })
    }

    pub fn to_bytes(&self) -> Vec<u8> {
        // converting a function to bytes???
        let mut res = vec![];
        res.extend(self.name.clone().into_bytes());
        res
    }
}

pub struct Ops<F: PrimeField + From<i32>> {
    pub add: Op<F>,
    pub mul: Op<F>,
}

impl<F: PrimeField + From<i32>> Ops<F> {
    pub fn new() -> Ops<F> {
        Ops {
            add: Op::new(
                "add".to_string(),
                |a, b| a + b,
                |mles| {
                    mles.iter().fold(0, |mut init, mle| {
                        init = max(init, mle.degree());
                        init
                    })
                },
            ),
            mul: Op::new(
                "mul".to_string(),
                |a, b| a * b,
                |mles| {
                    mles.iter().fold(0, |mut init, mle| {
                        init += mle.degree();
                        init
                    })
                },
            ),
        }
    }
}
