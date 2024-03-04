use std::collections::HashMap;

use super::constraint_rep::Constraint;
use ark_ff::PrimeField;

#[derive(Debug, Clone)]
pub struct ConstraintSystem<F: PrimeField> {
    pub constraints: Vec<Constraint<F>>,
    pub total_constraint: usize,
    pub total_vars: usize,
    /// TODO: is this neccesary? since the sign can be over written
    /// returns the value(sign) of a variable
    /// the usize is the variable
    /// the F is the value(sign)
    // pub vars: BTreeMap<usize, F>,
    /// Returns the index of a value in the witness
    pub constant_map: HashMap<F, usize>,
}

impl<F: PrimeField> ConstraintSystem<F> {
    pub fn new(constraints: &Vec<Constraint<F>>) -> Self {
        Self {
            constraints: constraints.to_vec(),
            total_constraint: constraints.len(),
            total_vars: 0,
            constant_map: HashMap::new(),
        }
    }

    pub fn update_total_vars(&mut self) {
        let mut max_var = 0;

        for constraint in &self.constraints {
            for i in 0..constraint.a.len() {
                let (var, _) = constraint.a[i];
                if var > max_var {
                    max_var = var;
                };
            }

            for i in 0..constraint.b.len() {
                let (var, _) = constraint.b[i];
                if var > max_var {
                    max_var = var;
                };
            }

            for i in 0..constraint.c.len() {
                let (var, _) = constraint.c[i];
                if var > max_var {
                    max_var = var;
                };
            }
        }

        self.total_vars = max_var;
    }

    pub fn update_constant_map(&mut self) {
        let mut constants = self.total_vars + 1;

        self.constant_map.insert(F::one(), 0);
        self.constant_map.insert(F::one().neg(), constants);
        constants = constants + 1;
        self.constant_map.insert(F::zero(), constants);
        constants = constants + 1;

        for constraint in &self.constraints {
            for i in 0..constraint.a.len() {
                let (_, val) = constraint.a[i];
                match self.constant_map.get(&val) {
                    Some(_) => {}
                    None => {
                        self.constant_map.insert(val, constants);
                        constants = constants + 1;
                    }
                }
            }

            for i in 0..constraint.b.len() {
                let (_, val) = constraint.b[i];
                match self.constant_map.get(&val) {
                    Some(_) => {}
                    None => {
                        self.constant_map.insert(val, constants);
                        constants = constants + 1;
                    }
                }
            }

            for i in 0..constraint.c.len() {
                let (_, val) = constraint.c[i];
                match self.constant_map.get(&val) {
                    Some(_) => {}
                    None => {
                        self.constant_map.insert(val, constants);
                        constants = constants + 1;
                    }
                }
            }
        }
    }

    pub fn create_constraints(
        &self,
        a: Option<(usize, F)>,
        b: Option<(usize, F)>,
        c: Option<(usize, F)>,
    ) -> Constraint<F> {
        let a = match a {
            Some(val) => vec![val],
            None => vec![],
        };

        let b = match b {
            Some(val) => vec![val],
            None => vec![],
        };

        let c = match c {
            Some(val) => vec![val],
            None => vec![],
        };

        Constraint::new_without_op(&a, &b, &c)
    }

    pub fn simplify_constraint_system(&mut self) {
        let mut res = vec![];

        self.update_total_vars();

        let mut new_constraints_from_a = vec![];
        let mut new_constraints_from_b = vec![];
        let mut new_constraints_from_c = vec![];

        for i in 0..self.constraints.len() {
            self.constraints[i].move_to_fill();

            let mut new_constraint = self.constraints[i].clone();

            while new_constraint.a.len() > 1
                || new_constraint.b.len() > 1
                || new_constraint.c.len() > 1
            {

                loop {
                    if new_constraint.a.len() <= 1 {
                        break;
                    }

                    self.total_vars = self.total_vars + 1;

                    if new_constraint.a.len() >= 2 {
                        let a = new_constraint.a.pop();
                        let b = new_constraint.a.pop();
                        let c = Some((self.total_vars, F::one()));
                        let constraint = self.create_constraints(a, b, c);

                        new_constraints_from_a.push(constraint);

                        new_constraint.a.push(c.unwrap());
                    }
                }

                loop {
                    if new_constraint.b.len() <= 1 {
                        break;
                    }

                    self.total_vars = self.total_vars + 1;

                    if new_constraint.b.len() >= 2 {
                        let a = new_constraint.b.pop();
                        let b = new_constraint.b.pop();
                        let c = Some((self.total_vars, F::one()));
                        let constraint = self.create_constraints(a, b, c);

                        new_constraints_from_b.push(constraint);

                        new_constraint.b.push(c.unwrap());
                    }
                }

                loop {
                    if new_constraint.c.len() <= 1 {
                        break;
                    }

                    self.total_vars = self.total_vars + 1;

                    if new_constraint.c.len() >= 2 {
                        let a = new_constraint.c.pop();
                        let b = new_constraint.c.pop();
                        let c = Some((self.total_vars, F::one()));
                        let constraint = self.create_constraints(a, b, c);

                        new_constraints_from_c.push(constraint);

                        new_constraint.c.push(c.unwrap());
                    }
                }
            }
            res.push(new_constraint);
        }

        res.extend(new_constraints_from_a);
        res.extend(new_constraints_from_b);
        res.extend(new_constraints_from_c);

        self.total_constraint = res.len();
        self.constraints = res;
        self.update_total_vars();
        self.update_constant_map();
    }
}

#[cfg(test)]
mod test {

    use crate::{circom_gkr::circom_read::read_circom_constraint, gkr::{self, circuit::Circuit, protocol::GKR}, polynomials::multilinear_poly::MultilinearPolynomialTrait};
    use ark_bn254::Fr;
    use ark_ff::Zero;
    
    #[test]
    fn test_simplifying() {
        let mut constraint_system =
        read_circom_constraint::<Fr>("src/circom_gkr/circom/test_constraints.json".to_string());

        constraint_system.simplify_constraint_system();
        
        // dbg!(&constraint_system.total_constraint);
        // dbg!(&constraint_system.constraints);
        // dbg!(&constraint_system.total_vars);
        // dbg!(&constraint_system.constant_map);

        let circuit: Circuit = constraint_system.into();

        // dbg!(&circuit.layers[0].len());
        // dbg!(&circuit.layers[1].len());
        // dbg!(&circuit.layers[2].len());
        // dbg!(&circuit.layers[3].len());

        let circuit_input = vec![
            // constants
            Fr::from(1),
            // intermediate vars
            Fr::from(5),
            Fr::from(3),
            Fr::from(2),
            Fr::from(1),
            Fr::from(1),
            Fr::from(5),
            Fr::from(3),
            Fr::from(10),
            Fr::from(7),
            // constants
            Fr::from(-1),
            Fr::from(0),
        ];

        // let circuit_input = vec![
        //     Fr::from(1),
        //     Fr::from(15),
        //     Fr::from(5),
        //     Fr::from(3),
        //     Fr::from(-1),
        //     Fr::from(0),  
        // ];


        let circuit_eval = circuit.evaluate(circuit_input);
       dbg!(&circuit.w_mle(circuit_eval[0].clone()).number_of_vars());
       dbg!(&circuit.w_mle(circuit_eval[1].clone()).number_of_vars());
       dbg!(&circuit.w_mle(circuit_eval[2].clone()).number_of_vars());
       dbg!(&circuit.w_mle(circuit_eval[3].clone()).number_of_vars());

        // let gkr_proof = GKR::prove(circuit.clone(), circuit_input.clone());
        // let gkr_verify = GKR::verify(circuit_input, gkr_proof.clone(), circuit).unwrap();
        // assert!(gkr_verify);
        // assert!(gkr_proof.sumcheck_proofs[0].sum == Fr::zero(), "Invalid circuit eval");

    }
}
