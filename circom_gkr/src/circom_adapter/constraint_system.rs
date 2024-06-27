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

    // Keeps track of the total number of variables in the constraint system
    // by keeping track of the highest variable addressed
    pub fn update_total_vars(&mut self) {
        let mut max_var = 0;

        for constraint in &self.constraints {
            for i in 0..constraint.a.len() {
                let (var, _) = constraint.a[i];
                if var > max_var {
                    max_var = var;
                    dbg!("{}", max_var);
                };
            }
            
            for i in 0..constraint.b.len() {
                let (var, _) = constraint.b[i];
                if var > max_var {
                    max_var = var;
                    dbg!("{}", max_var);
                };
            }
            
            for i in 0..constraint.c.len() {
                let (var, _) = constraint.c[i];
                if var > max_var {
                    max_var = var;
                    dbg!("{}", max_var);
                };
            }
        }

        self.total_vars = max_var;
    }

    // Keeps track of the new variables created
    // These new variables are added before -1 and 0
    pub fn update_constant_map(&mut self) {

        for constraint in &self.constraints {
            for i in 0..constraint.a.len() {
                let (_, val) = constraint.a[i];

                if val == F::one() || val == F::one().neg() {
                    continue;
                }

                match self.constant_map.get(&val) {
                    Some(_) => {}
                    None => {
                        self.total_vars = self.total_vars + 1;

                        println!("Inserting {}: {val}", self.total_vars);
                        
                        self.constant_map.insert(val, self.total_vars);
                    }
                }
            }
            
            for i in 0..constraint.b.len() {
                let (_, val) = constraint.b[i];

                if val == F::one() || val == F::one().neg() {
                    continue;
                }

                match self.constant_map.get(&val) {
                    Some(_) => {}
                    None => {
                        self.total_vars = self.total_vars + 1;

                        println!("Inserting {}: {val}", self.total_vars);
                        
                        self.constant_map.insert(val, self.total_vars);
                    }
                }
            }
            
            for i in 0..constraint.c.len() {
                let (_, val) = constraint.c[i];

                if val == F::one() || val == F::one().neg() {
                    continue;
                }

                match self.constant_map.get(&val) {
                    Some(_) => {}
                    None => {
                        self.total_vars = self.total_vars + 1;

                        println!("Inserting {}: {val}", self.total_vars);

                        self.constant_map.insert(val, self.total_vars);
                    }
                }
            }
        }

        println!("Total vars before constants {}", self.total_vars);
        
    }

    pub fn finalize_constants_generation(&mut self) {
        self.update_constant_map();

        println!("{:?}", self.constant_map);
        println!("{:?}", self.constant_map.len());

        self.constant_map.insert(F::one(), 0);

        match self.constant_map.get(&F::one().neg()) {
            Some(_) => {}
            None => {
                self.total_vars = self.total_vars + 1;
                self.constant_map.insert(F::one().neg(), self.total_vars);
            }
        }

        match self.constant_map.get(&F::zero()) {
            Some(_) => {}
            None => {
                self.total_vars = self.total_vars + 1;
                self.constant_map.insert(F::zero(), self.total_vars);
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
        // self.update_total_vars();
        // self.update_constant_map();
        self.finalize_constants_generation();
    }
}

#[cfg(test)]
mod test {

    use crate::circom_adapter::circom_read::read_circom_constraint;
    use ark_bn254::Fr;
    use ark_ff::Zero;
    use gkr::{circuit::Circuit, protocol::GKR};

    #[test]
    fn test_simplifying() {
        let mut constraint_system =
            read_circom_constraint::<Fr>("other_assets/test_constraints.json".to_string());

        constraint_system.simplify_constraint_system();

        let circuit: Circuit = constraint_system.into();

        let circuit_input = vec![
            // constants
            Fr::from(1), // 0
            // intermediate vars
            Fr::from(5),  // 1
            Fr::from(3),  // 2
            Fr::from(2),  // 3
            Fr::from(4),  // 4
            Fr::from(16), // 5
            Fr::from(80), // 6
            Fr::from(12), // 7
            Fr::from(94), // 8
            Fr::from(82), // 9
            Fr::from(0), // 10
            Fr::from(0), // 11
            // constants
            Fr::from(-1), // 12
            Fr::from(0),  // 13
            // new variable created
        ];


        // let circuit_eval = circuit.evaluate(circuit_input.clone());

        // dbg!(&circuit_eval);

        //    dbg!(&circuit.w_mle(circuit_eval[0].clone()).number_of_vars());
        //    dbg!(&circuit.w_mle(circuit_eval[1].clone()).number_of_vars());
        //    dbg!(&circuit.w_mle(circuit_eval[2].clone()).number_of_vars());
        //    dbg!(&circuit.w_mle(circuit_eval[3].clone()).number_of_vars());

        let gkr_proof = GKR::prove(&circuit, &circuit_input);

        let gkr_verify = GKR::verify(circuit_input, gkr_proof.clone(), circuit).unwrap();

        dbg!(gkr_verify);

        dbg!(&gkr_proof.sumcheck_proofs[0].sum);

        assert!(
            gkr_proof.sumcheck_proofs[0].sum == Fr::zero(),
            "Invalid circuit eval"
        );
    }

    #[test]
    fn test_generate_intermediate_variables() {
        let mut constraint_system =
            read_circom_constraint::<Fr>("other_assets/test_constraints.json".to_string());

            // Constraints generated
            // Constraint 1: -x * x = -w => (x * -1) * (x * 1) = (w * -1)
            // Constraint 2: -a * w = -y => (a * -1) * (w * 1) = (y * -1)
            // Constraint 3: -b * x = -z => (b * -1) * (x * 1) = (z * -1)
            // Constraint 4: c + y + z -n = 0 => (c * 1) + (y * 1) + (z * 1) + (n * -1)
            // New constraint is to be created from Constraint 4
            // Constraint 5: n - z = newvar
            // Constraint 6: y + c = newvar
        constraint_system.simplify_constraint_system();

        // The main 3 values in the map are 1, -1, 0
        // They are essential for the circuit formation
        assert!(constraint_system.constant_map.len() == 3, "Wrong number of constants");

        // Check if the constants are placed at their right index

    }

}
