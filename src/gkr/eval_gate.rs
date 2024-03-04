use std::{fmt::Debug, ops::Add};

use ark_ff::{BigInteger, PrimeField};
use ark_std::iterable::Iterable;

use crate::{
    polynomials::multilinear_poly::{MultilinearPolynomial, MultilinearPolynomialTrait},
    utils::get_binary_string,
};

/// Evalgate is eqivalent to f(b, c)
/// f(b, c) = add_i(a, b, c)(w_mle(b) + w_mle(c)) + add_i(a, b, c)(w_mle(b) * w_mle(c))
#[derive(Debug, Clone)]
pub struct EvalGate<F: PrimeField> {
    pub r: Vec<F>,
    pub add_i_mle: Vec<MultilinearPolynomial<F>>,
    pub mul_i_mle: Vec<MultilinearPolynomial<F>>,
    pub w_b_mle: Vec<MultilinearPolynomial<F>>,
    pub w_c_mle: Vec<MultilinearPolynomial<F>>,
}

impl<F: PrimeField> EvalGate<F> {
    pub fn new(
        r: &Vec<F>,
        add_i_mle: &Vec<MultilinearPolynomial<F>>,
        mul_i_mle: &Vec<MultilinearPolynomial<F>>,
        w_b_mle: &Vec<MultilinearPolynomial<F>>,
        w_c_mle: &Vec<MultilinearPolynomial<F>>,
    ) -> Self {
        Self {
            r: r.clone(),
            add_i_mle: add_i_mle.clone(),
            mul_i_mle: mul_i_mle.clone(),
            w_b_mle: w_b_mle.clone(),
            w_c_mle: w_c_mle.clone(),
        }
    }
}

impl<F: PrimeField> MultilinearPolynomialTrait<F> for EvalGate<F> {
    fn partial_eval(&self, x: &Vec<(usize, F)>) -> EvalGate<F> {

        let bnc = x.iter().fold(vec![vec![]; 2], |mut acc, (index, ele)| {
            if *index < self.w_b_num_of_vars() {
                acc[0].push((index.clone(), *ele));
                acc
            } else {
                acc[1].push((index.clone(), *ele));
                acc
            }
        });

        let b_eval_points = bnc[0].clone();
        let c_eval_points: Vec<(usize, F)> = bnc[1]
            .clone()
            .into_iter()
            .map(|(index, coeff)| (index - self.w_b_num_of_vars(), coeff))
            .collect();

        let eval_domain: Vec<(usize, F)> = x
            .clone()
            .into_iter()
            .map(|(index, point)| (self.r.len() + index, point))
            .collect();

        let w_b_mle = self
            .w_b_mle
            .clone()
            .into_iter()
            .fold(vec![], |mut acc, poly| {
                acc.push(poly.partial_eval(&b_eval_points));
                acc
            });
        let w_c_mle = self
            .w_c_mle
            .clone()
            .into_iter()
            .fold(vec![], |mut acc, poly| {
                acc.push(poly.partial_eval(&c_eval_points));
                acc
            });

        let add_i_mle = self
            .add_i_mle
            .clone()
            .into_iter()
            .fold(vec![], |mut acc, poly| {
                acc.push(poly.partial_eval(&eval_domain).relabel());
                acc
            });

        let mul_i_mle = self
            .mul_i_mle
            .clone()
            .into_iter()
            .fold(vec![], |mut acc, poly| {
                acc.push(poly.partial_eval(&eval_domain).relabel());
                acc
            });

        EvalGate::new(&self.r.clone(), &add_i_mle, &mul_i_mle, &w_b_mle, &w_c_mle)
    }

    fn evaluate(&self, x: &Vec<(usize, F)>) -> F {
        let mut eval_domain: Vec<(usize, F)> = vec![];

        let (b, c) = x.split_at(self.w_b_num_of_vars());

        eval_domain = self.r.clone().into_iter().enumerate().collect();
        eval_domain.extend(b);
        eval_domain.extend(c);

        eval_domain = eval_domain
            .into_iter()
            .fold(vec![], |mut acc, (_, ele)| {
                acc.push(ele);
                acc
            })
            .into_iter()
            .enumerate()
            .collect();

        let b_eval_points = b
            .iter()
            .fold(vec![], |mut acc, (_, ele)| {
                acc.push(*ele);
                acc
            })
            .into_iter()
            .enumerate()
            .collect();

        let c_eval_points = c
            .iter()
            .fold(vec![], |mut acc, (_, ele)| {
                acc.push(*ele);
                acc
            })
            .into_iter()
            .enumerate()
            .collect();

        let mut res = F::zero();
            // dbg!(&eval_domain);
            // dbg!(&self.mul_i_mle[0].number_of_vars());
        for i in 0..self.w_b_mle.len() {

            let b_eval = self.w_b_mle[i].evaluate(&b_eval_points);
            let c_eval = self.w_c_mle[i].evaluate(&c_eval_points);

            let add_result = self.add_i_mle[i].evaluate(&eval_domain) * (b_eval + c_eval);
            let mul_result = self.mul_i_mle[i].evaluate(&eval_domain) * (b_eval * c_eval);

            res += add_result + mul_result;
        }

        res
    }

    fn number_of_vars(&self) -> usize {
        self.w_b_num_of_vars() + self.w_c_num_of_vars()        
    }

    fn to_bytes(&self) -> Vec<u8> {
        let res = self.r.iter().fold(vec![], |mut acc, field| {
            acc.extend(field.into_bigint().to_bytes_be());
            acc
        });

        let res = self.add_i_mle.iter().fold(res, |mut res, poly| {
            res.extend(poly.to_bytes());
            res
        });

        let res = self.mul_i_mle.iter().fold(res, |mut res, poly| {
            res.extend(poly.to_bytes());
            res
        });

        let res = self.w_b_mle.iter().fold(res, |mut res, poly| {
            res.extend(poly.to_bytes());
            res
        });

        let res = self.w_c_mle.iter().fold(res, |mut res, poly| {
            res.extend(poly.to_bytes());
            res
        });

        res
    }

    fn relabel(&self) -> Self {
        let new_add_i_mle = self
            .add_i_mle
            .clone()
            .into_iter()
            .fold(vec![], |mut acc, poly| {
                acc.push(poly.relabel());
                acc
            });

        let new_mul_i_mle = self
            .mul_i_mle
            .clone()
            .into_iter()
            .fold(vec![], |mut acc, poly| {
                acc.push(poly.relabel());
                acc
            });

        let new_w_b_mle = self
            .w_b_mle
            .clone()
            .into_iter()
            .fold(vec![], |mut acc, poly| {
                acc.push(poly.relabel());
                acc
            });

        let new_w_c_mle = self
            .w_c_mle
            .clone()
            .into_iter()
            .fold(vec![], |mut acc, poly| {
                acc.push(poly.relabel());
                acc
            });

        Self {
            r: self.r.clone(),
            add_i_mle: new_add_i_mle,
            mul_i_mle: new_mul_i_mle,
            w_b_mle: new_w_b_mle,
            w_c_mle: new_w_c_mle,
        }
    }

    fn additive_identity() -> Self {
        Self {
            r: vec![],
            add_i_mle: vec![],
            mul_i_mle: vec![],
            w_b_mle: vec![],
            w_c_mle: vec![],
        }
    }

    // Evaluates the sum over the boolean hypercube and returns the sum
    fn sum_over_the_boolean_hypercube(&self) -> F {
        let mut res = F::zero();
        if self.number_of_vars() == 0 {
            return F::zero();
        };
        let vars_len = self.number_of_vars();
        for i in 0..2_usize.pow(vars_len as u32) {
            let boolean_vec: Vec<F> = get_binary_string(i, vars_len)
                .chars()
                .into_iter()
                .map(|var| if var == '0' { F::zero() } else { F::one() })
                .collect();
            let eval_domain = boolean_vec.into_iter().enumerate().collect();
            res += self.evaluate(&eval_domain);
        }
        res
    }
}

impl<F: PrimeField> EvalGate<F> {
    fn w_b_num_of_vars(&self) -> usize {
        self.w_b_mle
            .first()
            .map(|val| val.number_of_vars())
            .unwrap_or(0)
    }

    fn w_c_num_of_vars(&self) -> usize {
        self.w_c_mle
            .first()
            .map(|val| val.number_of_vars())
            .unwrap_or(0)
    }
}

// Implement native addition for Evalgate
impl<F: PrimeField> Add for EvalGate<F> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self {
        let r = if self.r.is_empty() { rhs.r } else { self.r };

        let mut new_add_i_mle = self.add_i_mle.clone();
        new_add_i_mle.extend(rhs.add_i_mle);

        let mut new_mul_i_mle = self.mul_i_mle.clone();
        new_mul_i_mle.extend(rhs.mul_i_mle);

        let mut new_w_b_mle = self.w_b_mle.clone();
        new_w_b_mle.extend(rhs.w_b_mle);

        let mut new_w_c_mle = self.w_c_mle.clone();
        new_w_c_mle.extend(rhs.w_c_mle);

        Self {
            r,
            add_i_mle: new_add_i_mle,
            mul_i_mle: new_mul_i_mle,
            w_b_mle: new_w_b_mle,
            w_c_mle: new_w_c_mle,
        }
    }
}

#[cfg(test)]
mod test {
    use ark_ff::{Fp64, MontBackend, MontConfig};

    use crate::{
        gkr::circuit::{Circuit, Layer, Wire},
        polynomials::multilinear_poly::{MultilinearPolynomial, MultilinearPolynomialTrait},
    };

    use super::EvalGate;

    #[derive(MontConfig)]
    #[modulus = "17"]
    #[generator = "3"]
    pub struct FqConfig;
    pub type Fq = Fp64<MontBackend<FqConfig, 1>>;

    fn create_circuit() -> Circuit {
        let layer0 = Layer::new(vec![], vec![Wire::new(0, 0, 1)]);
        let layer1 = Layer::new(vec![Wire::new(0, 0, 1)], vec![Wire::new(1, 1, 2)]);
        let layer2 = Layer::new(
            vec![Wire::new(0, 0, 1), Wire::new(2, 4, 5)],
            vec![Wire::new(1, 2, 3)],
        );

        let new_circuit = Circuit::new(vec![layer0, layer1, layer2]);

        new_circuit
    }

    #[test]
    fn test_eval_gate_eval() {
        let circuit = create_circuit();

        let circuit_eval = circuit.evaluate(vec![
            Fq::from(5),
            Fq::from(2),
            Fq::from(3),
            Fq::from(4),
            Fq::from(9),
            Fq::from(8),
        ]);

        dbg!(circuit.depth);

        let layer_2_mle = circuit.layer_mle(2);

        // Test mul gate
        let r = vec![Fq::from(0), Fq::from(1)];
        let w_mle = MultilinearPolynomial::<Fq>::interpolate(circuit_eval[3].clone());

        let eval_gate = EvalGate::new(
            &r,
            &vec![layer_2_mle[0].clone()],
            &vec![layer_2_mle[1].clone()],
            &vec![w_mle.clone()],
            &vec![w_mle],
        );

        let b = vec![Fq::from(0), Fq::from(1), Fq::from(0)];
        let c = vec![Fq::from(0), Fq::from(1), Fq::from(1)];

        let mut eval_points = b;
        eval_points.extend(c);

        let eval_tup = eval_points
            .into_iter()
            .enumerate()
            .map(|(index, point)| (index + r.len(), point))
            .collect();

        let res = eval_gate.evaluate(&eval_tup);
        assert!(res == Fq::from(12));

        // Test add gate
        let r = vec![Fq::from(0), Fq::from(0)];

        let w_mle = MultilinearPolynomial::<Fq>::interpolate(circuit_eval[3].clone());
        let eval_gate = EvalGate::new(
            &r,
            &vec![layer_2_mle[0].clone()],
            &vec![layer_2_mle[1].clone()],
            &vec![w_mle.clone()],
            &vec![w_mle],
        );

        let b = vec![Fq::from(0), Fq::from(0), Fq::from(0)];
        let c = vec![Fq::from(0), Fq::from(0), Fq::from(1)];

        let mut eval_points = b;
        eval_points.extend(c);

        let eval_tup = eval_points
            .into_iter()
            .enumerate()
            .map(|(index, point)| (index + r.len(), point))
            .collect();

        let res = eval_gate.evaluate(&eval_tup);
        assert!(res == Fq::from(7));
    }

    #[test]
    fn test_eval_gate_partial_eval() {
        let circuit = create_circuit();

        let circuit_eval = circuit.evaluate(vec![
            Fq::from(5),
            Fq::from(2),
            Fq::from(3),
            Fq::from(4),
            Fq::from(9),
            Fq::from(8),
        ]);

        let [add_i_mle, mul_i_mle] = circuit.layer_mle(2);

        let r = vec![Fq::from(0), Fq::from(1)];
        let w_mle = MultilinearPolynomial::<Fq>::interpolate(circuit_eval[3].clone());

        let eval_gate = EvalGate::new(
            &r,
            &vec![add_i_mle],
            &vec![mul_i_mle],
            &vec![w_mle.clone()],
            &vec![w_mle.clone()],
        );

        let b = vec![Fq::from(0), Fq::from(1), Fq::from(0)];
        let c = vec![Fq::from(0), Fq::from(1), Fq::from(1)];

        let mut eval_points = b.clone();
        eval_points.extend(c);

        let eval_tup = eval_points
            .into_iter()
            .enumerate()
            .map(|(index, val)| (index, val))
            .collect();

        let partial_res = eval_gate.partial_eval(&eval_tup).relabel();

        let res = partial_res.relabel().evaluate(&vec![]);

        assert!(res == Fq::from(12));
    }
}
