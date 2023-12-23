use ark_ff::{BigInteger, PrimeField};

use crate::polynomials::multilinear_poly::{MultilinearPolynomial, MultilinearPolynomialTrait};

#[derive(Debug, Clone)]
pub struct EvalGate<F: PrimeField, MPT: MultilinearPolynomialTrait<F>> {
    pub r: Vec<F>,
    pub add_i_mle: MultilinearPolynomial<F>,
    pub mul_i_mle: MultilinearPolynomial<F>,
    pub w_b_mle: MPT,
    pub w_c_mle: MPT,
}

impl<F: PrimeField, MPT: MultilinearPolynomialTrait<F> + Clone> EvalGate<F, MPT> {
    pub fn new(
        r: Vec<F>,
        add_i_mle: MultilinearPolynomial<F>,
        mul_i_mle: MultilinearPolynomial<F>,
        w_b_mle: MPT,
        w_c_mle: MPT,
    ) -> Self {
        Self {
            r,
            add_i_mle,
            mul_i_mle,
            w_b_mle,
            w_c_mle,
        }
    }
}

impl<F: PrimeField, MPT: MultilinearPolynomialTrait<F> + Clone> MultilinearPolynomialTrait<F>
    for EvalGate<F, MPT>
{
    fn partial_eval(&self, x: Vec<(usize, F)>) -> EvalGate<F, MPT> {
        let bnc = x.iter().fold(vec![vec![]; 2], |mut acc, (index, ele)| {
            if *index > self.r.len() && *index < self.r.len() + self.w_b_mle.number_of_vars() {
                acc[0].push(*ele);
                acc
            } else if *index > self.r.len() + self.w_b_mle.number_of_vars() {
                acc[1].push(*ele);
                acc
            } else {
                acc
            }
        });

        let b_eval_points: Vec<(usize, F)> = if bnc[0].len() < self.w_b_mle.number_of_vars() {
            let mut b = vec![F::zero(); self.w_b_mle.number_of_vars() - bnc[0].len()];
            b.extend(&bnc[0]);
            b.into_iter().enumerate().collect()
        } else {
            bnc[0].clone().into_iter().enumerate().collect()
        };

        let c_eval_points: Vec<(usize, F)> = if bnc[1].len() < self.w_c_mle.number_of_vars() {
            let mut c = vec![F::zero(); self.w_c_mle.number_of_vars() - bnc[1].len()];
            c.extend(&bnc[1]);
            c.into_iter().enumerate().collect()
        } else {
            bnc[1].clone().into_iter().enumerate().collect()
        };

        let w_b_mle = self.w_b_mle.partial_eval(b_eval_points);
        let w_c_mle = self.w_c_mle.partial_eval(c_eval_points);

        let add_i_mle = self.add_i_mle.partial_eval(x.clone());
        let mul_i_mle = self.mul_i_mle.partial_eval(x);

        EvalGate::new(self.r.clone(), add_i_mle, mul_i_mle, w_b_mle, w_c_mle)
    }

    // f(b, c) = add_i(a, b, c)(w_mle(b) + w_mle(c)) + add_i(a, b, c)(w_mle(b) + w_mle(c))
    fn evaluate(&self, x: Vec<(usize, F)>) -> F {
        let mut b: Vec<(usize, F)> = vec![];
        let mut c: Vec<(usize, F)> = vec![];

        if x.len() > self.r.len() {
            b = x[self.r.len()..self.w_b_mle.number_of_vars() + self.r.len()].to_vec();
            c = x[self.w_b_mle.number_of_vars() + self.r.len()..].to_vec();
        }

        let mut eval_domain: Vec<(usize, F)> = self.r.clone().into_iter().enumerate().collect();
        eval_domain.extend(&b);
        eval_domain.extend(&c);

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

        let b_eval = self.w_b_mle.evaluate(b_eval_points);
        let c_eval = self.w_c_mle.evaluate(c_eval_points);

        let add_result = self.add_i_mle.evaluate(eval_domain.clone()) * (b_eval + c_eval);
        let mul_result = self.mul_i_mle.evaluate(eval_domain.clone()) * (b_eval * c_eval);
        add_result + mul_result
    }

    fn number_of_vars(&self) -> usize {
        self.w_b_mle.number_of_vars() + self.w_c_mle.number_of_vars()
    }

    fn to_bytes(&self) -> Vec<u8> {
        let mut res = self.r.iter().fold(vec![], |mut acc, field| {
            acc.extend(field.into_bigint().to_bytes_be());
            acc
        });
        res.extend(self.add_i_mle.to_bytes());
        res.extend(self.mul_i_mle.to_bytes());
        res.extend(self.w_b_mle.to_bytes());
        res.extend(self.w_c_mle.to_bytes());

        res
    }

    fn relabel(&self) -> Self {
        Self {
            r: self.r.clone(),
            add_i_mle: self.add_i_mle.relabel(),
            mul_i_mle: self.mul_i_mle.relabel(),
            w_b_mle: self.w_b_mle.relabel(),
            w_c_mle: self.w_c_mle.relabel(),
        }
    }

    fn additive_identity() -> Self {
        Self {
            r: vec![],
            add_i_mle: MultilinearPolynomial::new(vec![]),
            mul_i_mle: MultilinearPolynomial::new(vec![]),
            w_b_mle: MPT::additive_identity(),
            w_c_mle: MPT::additive_identity(),
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

    pub fn create_circuit() -> Circuit {
        let mut new_circuit = Circuit::new(vec![], 0);
        let layer0 = Layer::new(
            vec![Wire::new(0, 0, 1), Wire::new(2, 4, 5)],
            vec![Wire::new(1, 2, 3)],
        );
        let layer1 = Layer::new(vec![Wire::new(0, 0, 1)], vec![Wire::new(1, 1, 2)]);
        let layer2 = Layer::new(vec![], vec![Wire::new(0, 0, 1)]);

        new_circuit.layers.push(layer0);
        new_circuit.layers.push(layer1);
        new_circuit.layers.push(layer2);

        new_circuit.depth = new_circuit.layers.len() + 1;

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

        let layer_0_mle = circuit.layer_mle(0);

        let r = vec![Fq::from(0), Fq::from(0), Fq::from(1)];
        let w_mle = MultilinearPolynomial::<Fq>::interpolate(circuit_eval[0].clone());

        let eval_gate = EvalGate::new(
            r.clone(),
            layer_0_mle[0].clone(),
            layer_0_mle[1].clone(),
            w_mle.clone(),
            w_mle,
        );

        let b = vec![Fq::from(0), Fq::from(1), Fq::from(0)];
        let c = vec![Fq::from(0), Fq::from(1), Fq::from(1)];

        let mut eval_points = r.clone();
        eval_points.extend(b);
        eval_points.extend(c);
        let eval_tup = eval_points.into_iter().enumerate().collect();

        let res = eval_gate.evaluate(eval_tup);
        assert!(res == Fq::from(12));

        // TODO: Test add gate
        // let r = vec![Fq::from(0), Fq::from(0), Fq::from(1)];
        // let w_mle = MultilinearPolynomial::<Fq>::interpolate(circuit_eval[0].clone());

        // let eval_gate = EvalGate::new(r, layer_0_mle[0].clone(), layer_0_mle[1].clone(), w_mle);

        // let b = vec![Fq::from(0), Fq::from(1), Fq::from(0)];
        // let c = vec![Fq::from(0), Fq::from(1), Fq::from(1)];

        // let mut eval_points = b;
        // eval_points.extend(c);
        // let eval_tup = eval_points.into_iter().enumerate().collect();

        // let res = eval_gate.evaluate(eval_tup);
        // dbg!(&res);
        // assert!(res == Fq::from(12));
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

        let layer_0_mle = circuit.layer_mle(0);

        let r = vec![Fq::from(0), Fq::from(0), Fq::from(1)];
        let w_mle = MultilinearPolynomial::<Fq>::interpolate(circuit_eval[0].clone());

        let eval_gate = EvalGate::new(
            r.clone(),
            layer_0_mle[0].clone(),
            layer_0_mle[1].clone(),
            w_mle.clone(),
            w_mle,
        );

        let b = vec![Fq::from(0), Fq::from(1), Fq::from(0)];
        let c = vec![Fq::from(0), Fq::from(1), Fq::from(1)];

        let mut eval_points = b.clone();
        eval_points.extend(c);

        let eval_tup = eval_points
            .into_iter()
            .enumerate()
            .map(|(index, val)| (eval_gate.r.len() + index, val))
            .collect();

        let partial_res = eval_gate.partial_eval(eval_tup);
        let res = partial_res
            .relabel()
            .evaluate(r.into_iter().enumerate().collect());
        assert!(res == Fq::from(12));
    }
}
