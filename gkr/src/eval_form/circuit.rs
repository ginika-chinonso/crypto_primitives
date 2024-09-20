use ark_ff::PrimeField;
use ark_std::log2;
use polynomials::multilinear_polynomial::{
    eval_form::MLE,
    universal_mle::{ops::Op, universal_mle::UniversalMLE},
};

pub struct Circuit<F: PrimeField> {
    pub depth: usize,
    pub layers: Vec<Layer<F>>,
}

pub struct Gate<F: PrimeField> {
    pub output: usize,
    pub left: usize,
    pub right: usize,
    pub op: Op<F>,
}

impl<F: PrimeField> Gate<F> {
    pub fn new(output: usize, left: usize, right: usize, op: &Op<F>) -> Self {
        Self {
            output,
            left,
            right,
            op: op.clone(),
        }
    }

    // The inputs is meant to be a vector of field elements passed from the layer beneath
    // All gates in a layer share the same inputs
    pub fn evaluate(&self, inputs: &Vec<F>, outputs: &mut Vec<F>) {
        outputs[self.output] = (self.op.op)(inputs[self.left], inputs[self.right])
    }
}

pub struct Layer<F: PrimeField> {
    pub gates: Vec<Gate<F>>,
    // pub num_of_vars: usize,
    pub layer_ops: Vec<Op<F>>,
}

impl<F: PrimeField> Layer<F> {
    pub fn new(gates: Vec<Gate<F>>) -> Self {
        let mut layer_ops = vec![];

        for gate in &gates {
            if layer_ops.contains(&gate.op) {
                continue;
            } else {
                layer_ops.push(gate.op.clone());
            }
        }

        Self { gates, layer_ops }
    }

    pub fn num_of_vars(&self) -> usize {
        self.layer_output_num_of_vars() + (2 * self.layer_input_num_of_vars())
    }

    pub fn layer_output_num_of_vars(&self) -> usize {
        if self.gates.len() == 1 {
            1
        } else {
            log2(self.gates.len().next_power_of_two()) as usize
        }
    }

    pub fn layer_input_num_of_vars(&self) -> usize {
        (log2(self.gates.len().next_power_of_two()) + 1) as usize
    }
}

impl<F: PrimeField> Circuit<F> {
    // The layer at index 0 is the output layer
    // while the layer at index depth is the input layer
    pub fn new(layers: Vec<Layer<F>>) -> Self {
        Self {
            depth: layers.len(),
            layers,
        }
    }

    pub fn evaluate(&self, inputs: &Vec<F>) -> Vec<Vec<F>> {
        // ensure that the input length is correct

        let mut res = vec![vec![]; self.depth + 1];

        res[self.depth] = inputs.to_vec();

        let mut depth = self.depth;

        for _ in 0..self.layers.len() {
            let mut layer_evaluation = vec![F::zero(); self.layers[depth - 1].gates.len()];

            for j in 0..self.layers[depth - 1].gates.len() {
                self.layers[depth - 1].gates[j].evaluate(&res[depth], &mut layer_evaluation)
            }

            depth -= 1;

            res[depth] = layer_evaluation;
        }

        res
    }

    pub fn get_layer_op_mle(&self, layer: usize, op: Op<F>) -> UniversalMLE<F> {
        let layer = &self.layers[layer];

        let mut val = vec![F::zero(); 2_usize.pow(layer.num_of_vars() as u32)];

        for gate in &layer.gates {
            if gate.op.name == op.name {
                let output = format!(
                    "{:0width$b}",
                    gate.output,
                    width = layer.layer_output_num_of_vars()
                );
                let left = format!(
                    "{:0width$b}",
                    gate.left,
                    width = layer.layer_input_num_of_vars()
                );
                let right = format!(
                    "{:0width$b}",
                    gate.right,
                    width = layer.layer_input_num_of_vars()
                );
                let gate_binary_string = String::new() + &output + &left + &right;
                let index: usize = usize::from_str_radix(gate_binary_string.as_str(), 2).unwrap();
                val[index] = F::one();
            }
        }

        UniversalMLE::MLES(vec![UniversalMLE::Value(MLE::new(&val))], op)
    }

    pub fn get_layer_ops_mles(&self, layer: usize) -> Vec<UniversalMLE<F>> {
        self.layers[layer]
            .layer_ops
            .iter()
            .map(|op| self.get_layer_op_mle(layer, op.clone()))
            .collect()
    }
}

pub fn get_wmle<F: PrimeField>(evaluations: &Vec<F>) -> UniversalMLE<F> {
    let val = if evaluations.len() == 1 {
        let mut v = evaluations.clone();
        v.push(F::zero());
        v
    } else if evaluations.len().is_power_of_two() {
        evaluations.clone()
    } else {
        let mut v = evaluations.clone();
        v.append(&mut vec![F::zero(); v.len().next_power_of_two() - v.len()]);
        v
    };
    UniversalMLE::Value(MLE::new(&val))
}

#[cfg(test)]
pub mod tests {
    use crate::eval_form::circuit::Layer;

    use super::{get_wmle, Circuit, Gate};
    use ark_bn254::Fq;
    use ark_ff::PrimeField;
    use polynomials::multilinear_polynomial::{
        traits::MultilinearPolynomialTrait, universal_mle::ops::Ops,
    };

    pub fn create_circuit<F: PrimeField>() -> Circuit<F> {
        let ops = Ops::new();

        let layer_2 = Layer::new(vec![
            Gate::new(0, 0, 1, &ops.add),
            Gate::new(1, 2, 3, &ops.mul),
            Gate::new(2, 4, 5, &ops.add),
        ]);

        let layer_1 = Layer::new(vec![
            Gate::new(0, 0, 1, &ops.mul),
            Gate::new(1, 1, 2, &ops.mul),
        ]);

        let layer_0 = Layer::new(vec![Gate::new(0, 0, 1, &ops.add)]);

        Circuit::new(vec![layer_0, layer_1, layer_2])
    }

    #[test]
    pub fn test_create_gate() {
        let inputs = vec![Fq::from(5), Fq::from(7), Fq::from(9)];

        let mut outputs = vec![Fq::from(0); 3];

        let gate = Gate::<Fq>::new(1, 2, 1, &Ops::new().add);

        gate.evaluate(&inputs, &mut outputs);

        assert_eq!(outputs, vec![Fq::from(0), Fq::from(16), Fq::from(0)]);
    }

    #[test]
    pub fn test_create_circuit() {
        let inputs = vec![
            Fq::from(1),
            Fq::from(2),
            Fq::from(3),
            Fq::from(4),
            Fq::from(5),
            Fq::from(6),
        ];

        let circuit = create_circuit();

        let res = circuit.evaluate(&inputs);

        assert_eq!(
            res[3],
            vec![
                Fq::from(1),
                Fq::from(2),
                Fq::from(3),
                Fq::from(4),
                Fq::from(5),
                Fq::from(6),
            ],
            "Incorrect circuit evaluation"
        );
        assert_eq!(
            res[2],
            vec![Fq::from(3), Fq::from(12), Fq::from(11)],
            "Incorrect circuit evaluation"
        );
        assert_eq!(
            res[1],
            vec![Fq::from(36), Fq::from(132)],
            "Incorrect circuit evaluation"
        );
        assert_eq!(res[0], vec![Fq::from(168)], "Incorrect circuit evaluation");
    }

    #[test]
    pub fn test_get_layer_wmle() {
        let mut inputs = vec![
            Fq::from(1),
            Fq::from(2),
            Fq::from(3),
            Fq::from(4),
            Fq::from(5),
            Fq::from(6),
        ];

        inputs.extend(vec![
            Fq::from(0);
            inputs.len().next_power_of_two() - inputs.len()
        ]);

        let circuit = create_circuit();

        let circuit_evaluations = circuit.evaluate(&inputs);

        let layer_3_mle = get_wmle(&circuit_evaluations[3]);

        assert_eq!(
            layer_3_mle.evaluate(&vec![(1, Fq::from(0)), (2, Fq::from(0)), (3, Fq::from(0))]),
            Fq::from(1),
            "Incorrect layer mle result"
        );
        assert_eq!(
            layer_3_mle.evaluate(&vec![(1, Fq::from(0)), (2, Fq::from(0)), (3, Fq::from(1))]),
            Fq::from(2),
            "Incorrect layer mle result"
        );
        assert_eq!(
            layer_3_mle.evaluate(&vec![(1, Fq::from(0)), (2, Fq::from(1)), (3, Fq::from(0))]),
            Fq::from(3),
            "Incorrect layer mle result"
        );
        assert_eq!(
            layer_3_mle.evaluate(&vec![(1, Fq::from(0)), (2, Fq::from(1)), (3, Fq::from(1))]),
            Fq::from(4),
            "Incorrect layer mle result"
        );
        assert_eq!(
            layer_3_mle.evaluate(&vec![(1, Fq::from(1)), (2, Fq::from(0)), (3, Fq::from(0))]),
            Fq::from(5),
            "Incorrect layer mle result"
        );
        assert_eq!(
            layer_3_mle.evaluate(&vec![(1, Fq::from(1)), (2, Fq::from(0)), (3, Fq::from(1))]),
            Fq::from(6),
            "Incorrect layer mle result"
        );
    }

    #[test]
    pub fn test_get_layer_op_mle() {
        let circuit: Circuit<Fq> = create_circuit();

        let layer_2_mul_mle = circuit.get_layer_op_mle(2, Ops::new().mul);

        assert_eq!(
            layer_2_mul_mle.sum_over_the_boolean_hypercube(),
            Fq::from(1),
            "Wrong number of multiplication gates"
        );
        assert_eq!(
            layer_2_mul_mle.number_of_vars(),
            8,
            "Wrong number of variables for layer"
        );
        assert_eq!(
            layer_2_mul_mle.evaluate(&vec![
                (1, Fq::from(0)),
                (2, Fq::from(1)),
                (3, Fq::from(0)),
                (4, Fq::from(1)),
                (5, Fq::from(0)),
                (6, Fq::from(0)),
                (7, Fq::from(1)),
                (8, Fq::from(1))
            ]),
            Fq::from(1),
            "Should evaluate to one on valid mul gate"
        );

        let layer_2_add_mle = circuit.get_layer_op_mle(2, Ops::new().add);

        assert_eq!(
            layer_2_add_mle.sum_over_the_boolean_hypercube(),
            Fq::from(2),
            "Wrong number of addition gates"
        );
        assert_eq!(
            layer_2_add_mle.number_of_vars(),
            8,
            "Wrong number of variables for layer"
        );
        assert_eq!(
            layer_2_add_mle.evaluate(&vec![
                (1, Fq::from(0)),
                (2, Fq::from(0)),
                (3, Fq::from(0)),
                (4, Fq::from(0)),
                (5, Fq::from(0)),
                (6, Fq::from(0)),
                (7, Fq::from(0)),
                (8, Fq::from(1))
            ]),
            Fq::from(1),
            "Should evaluate to one on valid add gate"
        );
        assert_eq!(
            layer_2_add_mle.evaluate(&vec![
                (1, Fq::from(1)),
                (2, Fq::from(0)),
                (3, Fq::from(1)),
                (4, Fq::from(0)),
                (5, Fq::from(0)),
                (6, Fq::from(1)),
                (7, Fq::from(0)),
                (8, Fq::from(1))
            ]),
            Fq::from(1),
            "Should evaluate to one on valid add gate"
        );
    }

    #[test]
    pub fn test_get_layer_ops_mles() {
        let circuit: Circuit<Fq> = create_circuit();

        let layer_ops_mles = circuit.get_layer_ops_mles(2);

        assert_eq!(layer_ops_mles.len(), 2, "Incorrect number of layer ops mle");
    }
}
