use std::collections::BTreeMap;

use super::constraint_rep::Constraint;
use ark_ff::PrimeField;

#[derive(Debug, Clone)]
pub struct ConstraintSystem<F: PrimeField> {
    pub constraints: Vec<Constraint<F>>,
    pub total_constraint: usize,
    pub total_vars: usize,
    pub vars: BTreeMap<usize, F>,
}

impl<F: PrimeField> ConstraintSystem<F> {
    pub fn new(constraints: &Vec<Constraint<F>>) -> Self {
        Self {
            constraints: constraints.to_vec(),
            total_constraint: constraints.len(),
            total_vars: 0,
            vars: BTreeMap::new(),
        }
    }

    pub fn get_total_vars(&self) -> usize {
        self.total_vars
    }

    pub fn update_vars(&mut self) {
        for constraint in &self.constraints {
            for i in 0..constraint.a.len() {
                let (var, val) = constraint.a[i];
                self.vars.insert(var, val);
            }

            for i in 0..constraint.b.len() {
                let (var, val) = constraint.b[i];
                self.vars.insert(var, val);
            }

            for i in 0..constraint.c.len() {
                let (var, val) = constraint.c[i];
                self.vars.insert(var, val);
            }
        }

        self.total_vars = self.get_last_var();
    }

    pub fn get_last_var(&self) -> usize {
        match self.vars.keys().cloned().max() {
            Some(num) => return num,
            None => return 0,
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

        self.update_vars();

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
                dbg!(&new_constraint.a.len());
                dbg!(&new_constraint.b.len());
                dbg!(&new_constraint.c.len());

                loop {
                    if new_constraint.a.len() <= 1 {
                        break;
                    }

                    let i = self.get_last_var() + 1;

                    if new_constraint.a.len() >= 2 {
                        let a = new_constraint.a.pop();
                        let b = new_constraint.a.pop();
                        let c = Some((i, F::one()));
                        let constraint = self.create_constraints(a, b, c);

                        new_constraints_from_a.push(constraint);

                        new_constraint.a.push(c.unwrap());

                        self.vars.insert(i, F::one());

                    }
                }

                loop {
                    if new_constraint.b.len() <= 1 {
                        break;
                    }

                    let i = self.get_last_var() + 1;

                    if new_constraint.b.len() >= 2 {
                        let a = new_constraint.b.pop();
                        let b = new_constraint.b.pop();
                        let c = Some((i, F::one()));
                        let constraint = self.create_constraints(a, b, c);

                        dbg!(&constraint);
                        new_constraints_from_b.push(constraint);

                        new_constraint.b.push(c.unwrap());

                        self.vars.insert(i, F::one());

                    }

                }

                loop {
                    if new_constraint.c.len() <= 1 {
                        break;
                    }

                    let i = self.get_last_var() + 1;

                    if new_constraint.c.len() >= 2 {
                        let a = new_constraint.c.pop();
                        let b = new_constraint.c.pop();
                        let c = Some((i, F::one()));
                        let constraint = self.create_constraints(a, b, c);

                        new_constraints_from_c.push(constraint);

                        new_constraint.c.push(c.unwrap());

                        self.vars.insert(i, F::one());

                    }
                }
            }
            res.push(new_constraint);
        }

        res.extend(new_constraints_from_a);
        res.extend(new_constraints_from_b);
        res.extend(new_constraints_from_c);

        self.constraints = res.clone();
        self.total_constraint = res.len();
        self.update_vars();
    }
}

#[cfg(test)]
mod test {

    use crate::circom_gkr::circom_read::read_circom_constraint;
    use ark_bn254::Fr;


    #[test]
    fn test_simplifying() {
        let mut constraint_system =
            read_circom_constraint::<Fr>("src/circom_gkr/circom/test_constraints.json".to_string());

        constraint_system.simplify_constraint_system();

        dbg!(&constraint_system.constraints);
        dbg!(&constraint_system.total_constraint);
        dbg!(&constraint_system.total_vars);

        assert!(false);
    }
}
