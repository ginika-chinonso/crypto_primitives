use std::vec;

use ark_ff::PrimeField;

#[derive(Debug, Clone, PartialEq)]
pub struct Constraint<F: PrimeField> {
    pub a: Vec<(usize, F)>,
    pub b: Vec<(usize, F)>,
    pub c: Vec<(usize, F)>,
    pub op: Op,
}

#[derive(Debug, Clone, PartialEq)]
pub enum Op {
    ADD,
    MUL,
}

impl<F: PrimeField> Constraint<F> {
    pub fn new_with_op(
        a: &Vec<(usize, F)>,
        b: &Vec<(usize, F)>,
        c: &Vec<(usize, F)>,
        op: Op,
    ) -> Self {
        return Self {
            a: a.to_vec(),
            b: b.to_vec(),
            c: c.to_vec(),
            op,
        };
    }

    pub fn new_without_op(a: &Vec<(usize, F)>, b: &Vec<(usize, F)>, c: &Vec<(usize, F)>) -> Self {
        let mut a = a.clone();
        let mut b = b.clone();
        let mut c = c.clone();

        if a.is_empty() {
            if b.len() > 1 {
                let new_a = match b.pop() {
                    Some(val) => vec![val],
                    None => vec![],
                };
                a = new_a;
            }

            if c.len() > 1 {
                if a.is_empty() {
                    let new_a = match c.pop() {
                        Some((var, val)) => vec![(var, val.neg())],
                        None => vec![],
                    };

                    a = new_a;
                }
            };
        };

        if b.is_empty() {
            if a.len() > 1 {
                let new_b = match a.pop() {
                    Some(val) => vec![val],
                    None => vec![],
                };
                b = new_b;
            }

            if c.len() > 1 {
                if b.is_empty() {
                    let new_b = match c.pop() {
                        Some((var, val)) => vec![(var, val.neg())],
                        None => vec![],
                    };

                    b = new_b;
                }
            };
        };

        if c.is_empty() {
            if a.len() > 1 {
                let new_c = match a.pop() {
                    Some((var, val)) => vec![(var, val.neg())],
                    None => vec![],
                };
                c = new_c;
            }

            if b.len() > 1 {
                if a.is_empty() {
                    let new_c = match b.pop() {
                        Some((var, val)) => vec![(var, val.neg())],
                        None => vec![],
                    };

                    c = new_c;
                }
            };
        };

        return Self {
            a,
            b,
            c,
            op: Op::ADD,
        };
    }

    pub fn can_move(&self) -> bool {
        if self.a.is_empty() || self.b.is_empty() || self.c.is_empty() {
            if self.op == Op::MUL {
                return false;
            } else {
                return true;
            }
        } else {
            return false;
        }
    }

    pub fn is_anb_empty(&self) -> bool {
        if self.a.is_empty() && self.b.is_empty() {
            return true;
        } else {
            return false;
        }
    }

    pub fn len(&self) -> usize {
        self.a.len() + self.b.len() + self.c.len()
    }

    pub fn can_simplify(&self) -> bool {
        if self.len() > 3 {
            return true;
        } else {
            if self.a.len() > 1 || self.b.len() > 1 || self.c.len() > 1 {
                return true;
            } else {
                return false;
            }
        }
    }

    pub fn move_to_fill(&mut self) {
        if self.can_move() {
            let new_self = Constraint::new_without_op(&self.a, &self.b, &self.c);
            self.a = new_self.a;
            self.b = new_self.b;
            self.c = new_self.c;
            self.op = new_self.op;
        }
    }
}

#[cfg(test)]
mod test {

    use super::Op;

    use super::Constraint;
    use ark_bls12_381::Fr;

    pub fn create_simplified_constraint() -> Constraint<Fr> {
        let s = vec![
            Fr::from(1),
            Fr::from(3),
            Fr::from(2),
            Fr::from(9),
            Fr::from(20),
        ];

        let a = vec![(0_usize, s[0])];
        let b = vec![(4_usize, s[4])];
        let c = vec![(3_usize, s[3])];

        Constraint::new_with_op(&a, &b, &c, Op::MUL)
    }

    pub fn create_can_move_constraint() -> Constraint<Fr> {
        let s = vec![
            Fr::from(1),
            Fr::from(3),
            Fr::from(2),
            Fr::from(9),
            Fr::from(20),
        ];

        let a = vec![(0_usize, s[0]), (4_usize, s[4])];
        let b = vec![];
        let c = vec![(3_usize, s[3]), (3_usize, s[3])];

        Constraint::new_with_op(&a, &b, &c, Op::ADD)
    }

    #[test]
    pub fn test_can_simplify() {
        let constraint = create_simplified_constraint();

        let is_simplifiable = constraint.can_simplify();

        assert!(is_simplifiable == false, "Can simplify check failed");
    }

    #[test]
    pub fn test_constraint_vars_len() {
        let constraint = create_simplified_constraint();

        assert!(constraint.len() == 3, "Incorrect constraint length");
    }

    #[test]
    pub fn test_move_to_fill() {
        let mut constraint = create_can_move_constraint();

        assert!(constraint.can_move(), "Can move check failed");

        constraint.move_to_fill();

        let a = vec![(0_usize, Fr::from(1))];
        let b = vec![(4_usize, Fr::from(20))];
        let c = vec![(3_usize, Fr::from(9)), (3_usize, Fr::from(9))];

        let expected = Constraint {
            a,
            b,
            c,
            op: Op::ADD,
        };

        assert_eq!(constraint, expected, "Move to fill failed");
    }
}
