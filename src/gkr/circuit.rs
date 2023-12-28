// - Implement circuit representation
// - circuit evaluation
// - w mle’s
// - add_i and mult_i mle’s
// - TODO: Implement sum over the boolean hypercube for multilinear poly
// prevent the creation of an empty layer

use ark_ff::PrimeField;

use crate::{
    polynomials::multilinear_poly::{MultilinearMonomial, MultilinearPolynomial},
    utils::{check_bit, get_binary_string},
};

pub struct Circuit {
    pub layers: Vec<Layer>,
    pub depth: usize,
}

impl Circuit {
    pub fn new(layers: Vec<Layer>, depth: usize) -> Self {
        Self { layers, depth }
    }

    pub fn evaluate<F: PrimeField>(&self, inputs: Vec<F>) -> Vec<Vec<F>> {
        // TODO: inner vec capacity should be sum of no of add and mul gates
        let mut res = vec![vec![F::zero(); self.depth]; self.depth];
        res[0] = inputs;
        let mut depth = 1_usize;
        for i in self.layers.clone() {
            let mut new_layer = vec![F::zero(); i.add_gates.len() + i.mul_gates.len()];
            for Wire {
                output,
                left,
                right,
            } in &i.add_gates
            {
                assert!(*output <= res[depth - 1].len());
                assert!(*left <= res[depth - 1].len());
                assert!(*right <= res[depth - 1].len());
                new_layer[*output] = res[depth - 1][*left] + res[depth - 1][*right];
            }

            for Wire {
                output,
                left,
                right,
            } in &i.mul_gates
            {
                assert!(*output <= res[depth - 1].len());
                assert!(*left <= res[depth - 1].len());
                assert!(*right <= res[depth - 1].len());
                new_layer[*output] = res[depth - 1][*left] * res[depth - 1][*right];
            }
            res[depth] = new_layer;
            depth += 1;
        }
        res
    }

    pub fn w_mle<F: PrimeField>(&self, layer_eval: Vec<F>) -> MultilinearPolynomial<F> {
        MultilinearPolynomial::interpolate(layer_eval)
    }

    pub fn layer_mle<F: PrimeField>(&self, layer: usize) -> [MultilinearPolynomial<F>; 2] {
        let layer_len = self.layers[layer].add_gates.len() + self.layers[layer].mul_gates.len();
        let add_i_mle = add_i_mle(&self.layers[layer].add_gates, layer_len);
        let mul_i_mle = mul_i_mle(&self.layers[layer].mul_gates, layer_len);
        [add_i_mle, mul_i_mle]
    }


}


#[derive(Debug, Clone)]
pub struct Wire {
    pub output: usize,
    pub left: usize,
    pub right: usize,
}

impl Wire {
    pub fn new(output: usize, left: usize, right: usize) -> Self {
        Self {
            output,
            left,
            right,
        }
    }
}

#[derive(Debug, Clone)]
pub struct Layer {
    add_gates: Vec<Wire>,
    mul_gates: Vec<Wire>,
}

impl Layer {
    pub fn new(add_gates: Vec<Wire>, mul_gates: Vec<Wire>) -> Self {
        Self {
            add_gates,
            mul_gates,
        }
    }
}

/////////////////////////////////////////////////////////////////////////
// Get padded binary string of a decimal number
/////////////////////////////////////////////////////////////////////////
pub fn get_index_poly<F: PrimeField>(
    index: usize,
    max_bit_count: usize,
) -> MultilinearPolynomial<F> {
    let mut res = MultilinearPolynomial::new(vec![MultilinearMonomial::new(F::one(), vec![])]);
    let binary = get_binary_string(index, max_bit_count);
    for j in 0..binary.len() {
        let i_char = binary.chars().nth(j).unwrap();
        if i_char == '0' {
            let i_rep = check_bit(0);
            res = res * i_rep;
        };
        if i_char == '1' {
            let i_rep = check_bit(1);
            res = res * i_rep;
        }
    }
    res
}

pub fn selector_poly<F: PrimeField>(
    a: usize,
    b: usize,
    c: usize,
    layer_len: usize,
) -> MultilinearPolynomial<F> {
    let a_poly = get_index_poly(a, layer_len);
    let b_poly = get_index_poly(b, layer_len);
    let c_poly = get_index_poly(c, layer_len);

    a_poly * b_poly * c_poly
}

pub fn add_i_mle<F: PrimeField>(gates: &[Wire], layer_len: usize) -> MultilinearPolynomial<F> {
    let res = MultilinearPolynomial::new(vec![]);
    gates.iter().fold(res, |res, i| {
        res + selector_poly(i.output, i.left, i.right, layer_len)
    })
}

pub fn mul_i_mle<F: PrimeField>(gates: &[Wire], layer_len: usize) -> MultilinearPolynomial<F> {
    let res = MultilinearPolynomial::new(vec![]);
    gates.iter().fold(res, |res, i| {
        res + selector_poly(i.output, i.left, i.right, layer_len)
    })
}

#[cfg(test)]
mod test {
    use ark_ff::{Fp64, MontBackend, MontConfig};

    use crate::gkr::circuit::selector_poly;

    use super::{Circuit, Layer, Wire};

    #[derive(MontConfig)]
    #[modulus = "17"]
    #[generator = "3"]
    pub struct FqConfig;
    pub type Fq = Fp64<MontBackend<FqConfig, 1>>;

    fn create_circuit() -> Circuit {
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
    fn test_circuit_eval() {
        let new_circuit = create_circuit();

        let circuit_eval = new_circuit.evaluate(vec![
            Fq::from(5),
            Fq::from(2),
            Fq::from(3),
            Fq::from(4),
            Fq::from(9),
            Fq::from(8),
        ]);
        dbg!(&circuit_eval);
        assert!(
            circuit_eval
                == vec![
                    vec![
                        Fq::from(5),
                        Fq::from(2),
                        Fq::from(3),
                        Fq::from(4),
                        Fq::from(9),
                        Fq::from(8)
                    ],
                    vec![Fq::from(7), Fq::from(12), Fq::from(17)],
                    vec![Fq::from(19), Fq::from(204)],
                    vec![Fq::from(3876)]
                ],
            "Incorrect result"
        );
    }

    #[test]
    fn test_selector_poly() {
        let add_i_mle = selector_poly(1, 0, 1, 3);
        let res: Fq = add_i_mle.sum_over_the_boolean_hypercube();
        assert!(res == Fq::from(1), "Incorrect result");
    }

    #[test]
    fn test_selector_poly_at_random_inputs() {
        let a = 5_usize;
        let b = 8_usize;
        let c = 9_usize;
        let add_i_mle = selector_poly(a, b, c, 4);
        let res: Fq = add_i_mle.sum_over_the_boolean_hypercube();
        assert!(res == Fq::from(1), "Incorrect result");
    }

    #[test]
    fn test_mul_i_mle() {
        let circuit = create_circuit();
        let layer_0_mle = circuit.layer_mle::<Fq>(0);
        let layer_0_mul_gates_mle = layer_0_mle[1].clone();
        let res = layer_0_mul_gates_mle.sum_over_the_boolean_hypercube();
        assert!(res == Fq::from(1), "Wrong number of gates");
        let layer_1_mle = circuit.layer_mle::<Fq>(1);
        let layer_1_mul_gates_mle = layer_1_mle[1].clone();
        let res = layer_1_mul_gates_mle.sum_over_the_boolean_hypercube();
        assert!(res == Fq::from(1), "Wrong number of gates");
        let layer_2_mle = circuit.layer_mle::<Fq>(2);
        let layer_2_mul_gates_mle = layer_2_mle[1].clone();
        let res = layer_2_mul_gates_mle.sum_over_the_boolean_hypercube();
        assert!(res == Fq::from(1), "Wrong number of gates");
    }

    #[test]
    fn test_add_i_mle() {
        let circuit = create_circuit();
        let layer_0_mle = circuit.layer_mle::<Fq>(0);
        let layer_0_add_gates_mle = layer_0_mle[0].clone();
        let res = layer_0_add_gates_mle.sum_over_the_boolean_hypercube();
        assert!(res == Fq::from(2), "Wrong number of gates");
        let layer_1_mle = circuit.layer_mle::<Fq>(1);
        let layer_1_add_gates_mle = layer_1_mle[0].clone();
        let res = layer_1_add_gates_mle.sum_over_the_boolean_hypercube();
        assert!(res == Fq::from(1), "Wrong number of gates");
        let layer_2_mle = circuit.layer_mle::<Fq>(2);
        let layer_2_add_gates_mle = layer_2_mle[0].clone();
        let res = layer_2_add_gates_mle.sum_over_the_boolean_hypercube();
        assert!(res == Fq::from(0), "Wrong number of gates");
    }
}
