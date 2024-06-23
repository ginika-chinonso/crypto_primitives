// Constraint system to GKR circuit

// A constraint ==>
// pub struct Constraint<F: PrimeField> {
//     pub a: Vec<(usize, F)>,
//     pub b: Vec<(usize, F)>,
//     pub c: Vec<(usize, F)>,
//     pub op: Op,
// }

// A wire ==>
// pub struct Wire {
//     pub output: usize,
//     pub left: usize,
//     pub right: usize,
// }

// Generate a constraint system from a simplified circuit
// TODO: Should it be possible to generate a circuit from an unsimlified constraint system?

// Circuit rep per constraint

//                                    0
//                                    |
//                                   add
//                       /                      \
//                     c                         -c
//                     |                          |
//                     op                         mul
//                /          \              /           \
//              new_a       new_b        new_c           -1
//               |            |            |              |
//              mul          mul          mul            mul
//            /    \        /   \        /    \         /   \
//           a   a_val     b   b_val    c    c_val     1    -1

use super::constraint_rep::Op;
use super::constraint_system::ConstraintSystem;
use ark_ff::PrimeField;
use gkr::circuit::{Circuit, Layer, Wire};

impl<F: PrimeField> Into<Circuit> for ConstraintSystem<F> {
    fn into(self) -> Circuit {
        let mut circuit = Circuit::new(vec![Layer::new(vec![], vec![]); 3]);

        for i in 0..self.constraints.len() {
            // let mut i_circuit = Circuit::new(vec![]);
            // TODO: is there a better way to do this without cloning?
            let mut constraint = self.constraints[i].clone();

            assert!(
                constraint.a.len() <= 1 && constraint.b.len() <= 1 && constraint.c.len() <= 1,
                "Please simplify Constraint system"
            );

            // all gates on the input layer are mulgates
            // so no need for the add_gates vec
            let mut mul_gates = vec![];

            // Rationale for the None cases
            // if None is returned, it means the variable is empty
            // It is expected to be the additive or multiplicative inverse depending on the constraint operation (this is the case only for a or b)
            // In the case of c, it is always expected to be zero

            let (a, a_val_index) = match constraint.a.pop() {
                Some((var, val)) => {
                    let val_index = *self.constant_map.get(&val).unwrap();
                    (var, val_index)
                }
                None => match constraint.op {
                    Op::ADD => {
                        let zero_index = *self.constant_map.get(&F::zero()).unwrap();
                        (zero_index, zero_index)
                    }
                    Op::MUL => {
                        let one_index = *self.constant_map.get(&F::one()).unwrap();
                        (one_index, one_index)
                    }
                },
            };
            let (b, b_val_index) = match constraint.b.pop() {
                Some((var, val)) => {
                    let val_index = self.constant_map.get(&val).unwrap().clone();
                    (var, val_index)
                }
                None => match constraint.op {
                    Op::ADD => {
                        let zero_index = *self.constant_map.get(&F::zero()).unwrap();
                        (zero_index, zero_index)
                    }
                    Op::MUL => {
                        let one_index = *self.constant_map.get(&F::one()).unwrap();
                        (one_index, one_index)
                    }
                },
            };
            let (c, c_val_index) = match constraint.c.pop() {
                Some((var, val)) => {
                    let val_index = self.constant_map.get(&val).unwrap().clone();
                    (var, val_index)
                }

                None => {
                    let zero_index = *self.constant_map.get(&F::zero()).unwrap();
                    (zero_index, zero_index)
                }
            };

            let (one, minus_one) = {
                let one = *self.constant_map.get(&F::one()).unwrap();
                let minus_one = *self.constant_map.get(&F::one().neg()).unwrap();
                (one, minus_one)
            };

            let mut output_index = 4 * i;

            let a_wire = Wire::new(output_index, a, a_val_index);
            mul_gates.push(a_wire);
            output_index += 1;

            let b_wire = Wire::new(output_index, b, b_val_index);
            mul_gates.push(b_wire);
            output_index += 1;

            let c_wire = Wire::new(output_index, c, c_val_index);
            mul_gates.push(c_wire);
            output_index += 1;

            let ones_wire = Wire::new(output_index, one, minus_one);
            mul_gates.push(ones_wire);

            // input layer constructed
            // let input_layer = Layer::new(vec![], mul_gates);

            // TODO: is there a better way to do this?
            // Since circuit depth is always constant, two can be hardcoded
            circuit.layers[2].mul_gates.extend(mul_gates);

            // build next layer
            let mut add_gates = vec![];
            let mut mul_gates = vec![];

            let mut output_index = 2 * i;
            let mut input_index = 4 * i;

            let op_wire = Wire::new(output_index, input_index, input_index + 1);
            output_index += 1;
            input_index += 2;

            match constraint.op {
                Op::ADD => {
                    add_gates.push(op_wire);
                }
                Op::MUL => {
                    mul_gates.push(op_wire);
                }
            };

            let minus_c_wire = Wire::new(output_index, input_index, input_index + 1);
            mul_gates.push(minus_c_wire);

            circuit.layers[1].add_gates.extend(add_gates);
            circuit.layers[1].mul_gates.extend(mul_gates);

            // Build output layer
            // output gates are always going to be add gates
            // mulgates can be left out
            let input_index = 2 * i;

            let c_minus_c = Wire::new(i, input_index, input_index + 1);

            circuit.layers[0].add_gates.extend(vec![c_minus_c]);
        }

        circuit
    }
}
