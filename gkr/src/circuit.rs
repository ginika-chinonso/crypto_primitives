// - Implement circuit representation
// - circuit evaluation
// - w mle’s
// - add_i and mult_i mle’s
// - TODO: Implement sum over the boolean hypercube for multilinear poly
// prevent the creation of an empty layer

use ark_serialize::*;
use std::cmp::max;

use ark_ff::PrimeField;

use polynomials::multilinear_polynomial::coef_form::{MultilinearMonomial, MultilinearPolynomial};
use utils::{check_bit, get_binary_string};

#[derive(Debug, Clone, CanonicalSerialize)]
pub struct Circuit {
    pub layers: Vec<Layer>,
    pub depth: usize,
}

impl Circuit {
    pub fn new(layers: Vec<Layer>) -> Self {
        let depth = layers.len();
        Self { layers, depth }
    }

    pub fn depth(&self) -> usize {
        self.layers.len()
    }

    pub fn evaluate<F: PrimeField>(&self, inputs: &Vec<F>) -> Vec<Vec<F>> {
        let mut res = vec![vec![F::zero(); inputs.len()]; self.depth + 1];
        let mut depth = self.depth;
        res[depth] = inputs.to_vec();

        loop {
            if depth == 0 {
                return res;
            };

            let mut layer_eval = vec![
                F::zero();
                self.layers[depth - 1].add_gates.len()
                    + self.layers[depth - 1].mul_gates.len()
            ];

            for Wire {
                output,
                left,
                right,
            } in self.layers[depth - 1].add_gates.clone()
            {
                // assert!(output <= res[depth].len());
                // assert!(left <= res[depth].len());
                // assert!(right <= res[depth].len());

                layer_eval[output] = res[depth][left] + res[depth][right];
            }

            for Wire {
                output,
                left,
                right,
            } in self.layers[depth - 1].mul_gates.clone()
            {
                // assert!(output <= res[depth].len());
                // assert!(left <= res[depth].len());
                // assert!(right <= res[depth].len());

                layer_eval[output] = res[depth][left] * res[depth][right];
            }

            depth -= 1;
            res[depth] = layer_eval;
        }
    }

    pub fn w_mle<F: PrimeField>(&self, layer_eval: Vec<F>) -> MultilinearPolynomial<F> {
        MultilinearPolynomial::interpolate(&layer_eval)
    }

    pub fn layer_mle<F: PrimeField>(&self, layer: usize) -> [MultilinearPolynomial<F>; 2] {
        let add_i_mle = self.add_i_mle(layer);
        let mul_i_mle = self.mul_i_mle(layer);
        [add_i_mle, mul_i_mle]
    }

    pub fn add_i_mle<F: PrimeField>(&self, layer: usize) -> MultilinearPolynomial<F> {
        let res = MultilinearPolynomial::new(vec![]);
        let layer = &self.layers[layer];
        let gates = layer.add_gates.to_vec();
        let layer_len = layer.len();
        let layer_input_len = layer.max_input_index() + 1;
        gates.iter().fold(res, |res, i| {
            res + selector_poly(i.output, i.left, i.right, layer_len, layer_input_len)
        })
    }

    pub fn mul_i_mle<F: PrimeField>(&self, layer: usize) -> MultilinearPolynomial<F> {
        let res = MultilinearPolynomial::new(vec![]);
        let layer = &self.layers[layer];
        let gates = layer.mul_gates.to_vec();
        let layer_len = layer.len();
        let layer_input_len = layer.max_input_index() + 1;
        gates.iter().fold(res, |res, i| {
            res + selector_poly(i.output, i.left, i.right, layer_len, layer_input_len)
        })
    }

    pub fn to_bytes(&self) -> Vec<u8> {
        let mut res = vec![];

        res.append(&mut self.depth.to_be_bytes().to_vec());

        for layer in &self.layers {
            res.append(&mut layer.to_bytes());
        }

        res
    }
}

// impl Serialize for Circuit {
//     fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
//     where
//         S: serde::Serializer {

//         // let res = vec![];
//         let mut seq = serializer.serialize_seq(Some(self.depth))?;

//         // res.push(self.depth.to_be_bytes());

//         for layer in &self.layers {
//             // res.push(layer.serialize(serializer))
//             let _ = seq.serialize_element(&layer)?;
//         }

//         seq.end()
//     }
// }

#[derive(Debug, Clone, Copy, CanonicalSerialize)]
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

    pub fn to_bytes(&self) -> Vec<u8> {
        let mut res = vec![];

        res.append(&mut self.output.to_be_bytes().to_vec());
        res.append(&mut self.left.to_be_bytes().to_vec());
        res.append(&mut self.right.to_be_bytes().to_vec());

        res
    }
}

#[derive(Debug, Clone, CanonicalSerialize)]
pub struct Layer {
    pub add_gates: Vec<Wire>,
    pub mul_gates: Vec<Wire>,
}

impl Layer {
    pub fn new(add_gates: Vec<Wire>, mul_gates: Vec<Wire>) -> Self {
        Self {
            add_gates,
            mul_gates,
        }
    }

    pub fn len(&self) -> usize {
        self.add_gates.len() + self.mul_gates.len()
    }

    pub fn max_input_index(&self) -> usize {
        let mut max_index = 0;

        for wire in &self.add_gates {
            max_index = max(wire.left, max_index);
            max_index = max(wire.right, max_index);
        }
        for wire in &self.mul_gates {
            max_index = max(wire.left, max_index);
            max_index = max(wire.right, max_index);
        }
        max_index
    }

    pub fn to_bytes(&self) -> Vec<u8> {
        let mut res = vec![];

        for wire in &self.add_gates {
            res.append(&mut wire.to_bytes());
        }
        for wire in &self.mul_gates {
            res.append(&mut wire.to_bytes());
        }

        res
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

    let max_bit_count = format!("{:b}", max_bit_count.checked_sub(1).unwrap_or_default()).len();

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
    input_layer_len: usize,
) -> MultilinearPolynomial<F> {
    let a_poly = get_index_poly(a, layer_len);
    let b_poly = get_index_poly(b, input_layer_len);
    let c_poly = get_index_poly(c, input_layer_len);

    // dbg!(&a, &layer_len);
    // dbg!(&b, &input_layer_len);
    // dbg!(&c, &input_layer_len);
    // dbg!(&a_poly);
    // dbg!(&b_poly);
    // dbg!(&c_poly);
    a_poly * b_poly * c_poly
}

#[cfg(test)]
pub mod test {
    use ark_ff::{Fp64, MontBackend, MontConfig};

    use crate::circuit::selector_poly;

    use super::{Circuit, Layer, Wire};
    use polynomials::multilinear_polynomial::traits::MultilinearPolynomialTrait;

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
    pub fn test_circuit_eval() {
        let new_circuit = create_circuit();

        let circuit_eval = new_circuit.evaluate(&vec![
            Fq::from(5),
            Fq::from(2),
            Fq::from(3),
            Fq::from(4),
            Fq::from(9),
            Fq::from(8),
        ]);

        let expected_result = vec![
            vec![Fq::from(3876)],
            vec![Fq::from(19), Fq::from(204)],
            vec![Fq::from(7), Fq::from(12), Fq::from(17)],
            vec![
                Fq::from(5),
                Fq::from(2),
                Fq::from(3),
                Fq::from(4),
                Fq::from(9),
                Fq::from(8),
            ],
        ];

        assert!(circuit_eval == expected_result, "Incorrect result");
    }

    #[test]
    fn test_selector_poly() {
        let add_i_mle = selector_poly(1, 0, 1, 3, 6);
        let res: Fq = add_i_mle.sum_over_the_boolean_hypercube();
        assert!(res == Fq::from(1), "Incorrect result");
    }

    #[test]
    fn test_selector_poly_at_random_inputs() {
        let a = 5_usize;
        let b = 8_usize;
        let c = 9_usize;
        let add_i_mle = selector_poly(a, b, c, 5, 10);
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
        assert!(res == Fq::from(0), "Wrong number of gates");
        let layer_1_mle = circuit.layer_mle::<Fq>(1);
        let layer_1_add_gates_mle = layer_1_mle[0].clone();
        let res = layer_1_add_gates_mle.sum_over_the_boolean_hypercube();
        assert!(res == Fq::from(1), "Wrong number of gates");
        let layer_2_mle = circuit.layer_mle::<Fq>(2);
        let layer_2_add_gates_mle = layer_2_mle[0].clone();
        let res = layer_2_add_gates_mle.sum_over_the_boolean_hypercube();
        assert!(res == Fq::from(2), "Wrong number of gates");
    }
}
