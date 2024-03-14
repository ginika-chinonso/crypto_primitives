use super::{
    constraint_rep::{Constraint, Op},
    constraint_system::ConstraintSystem,
};
use ark_ff::PrimeField;
use ark_std::fmt::Debug;
use serde_json::Value;
use std::fs::File;
use std::io::Read;

pub fn read_circom_constraint<F: PrimeField + Debug>(path: String) -> ConstraintSystem<F> {
    let mut res = vec![];

    let mut file = File::open(path).expect("failed to open file");

    let mut content = String::new();

    file.read_to_string(&mut content)
        .expect("failed to read file");

    let parsed_result: Value = serde_json::from_str(content.as_str()).unwrap();

    let constraints = parsed_result["constraints"].as_array().unwrap().clone();

    for i in 0..constraints.len() {
        let a_circom_con = constraints[i].as_array().unwrap()[0].as_object().unwrap();
        let b_circom_con = constraints[i].as_array().unwrap()[1].as_object().unwrap();
        let c_circom_con = constraints[i].as_array().unwrap()[2].as_object().unwrap();

        let mut a = vec![];
        let mut b = vec![];
        let mut c = vec![];

        for (key, val) in a_circom_con {
            let new_key = key.parse::<usize>().unwrap();
            let new_val = match F::from_str(
                serde_json::from_value::<String>(val.clone())
                    .unwrap()
                    .as_str(),
            ) {
                Ok(val) => Ok(val),
                Err(_) => Err("Failed to convert number to field element"),
            };
            a.push((new_key, new_val.unwrap()));
        }

        for (key, val) in b_circom_con {
            let new_key = key.parse::<usize>().unwrap();
            let new_val = match F::from_str(
                serde_json::from_value::<String>(val.clone())
                    .unwrap()
                    .as_str(),
            ) {
                Ok(val) => Ok(val),
                Err(_) => Err("Failed to convert number to field element"),
            };
            b.push((new_key, new_val.unwrap()));
        }

        for (key, val) in c_circom_con {
            let new_key = key.parse::<usize>().unwrap();
            let new_val = match F::from_str(
                serde_json::from_value::<String>(val.clone())
                    .unwrap()
                    .as_str(),
            ) {
                Ok(val) => Ok(val),
                Err(_) => Err("Failed to convert number to field element"),
            };
            c.push((new_key, new_val.unwrap()));
        }

        let new_constraint = if a.is_empty() || b.is_empty() {
            Constraint {
                a,
                b,
                c,
                op: Op::ADD,
            }
        } else {
            Constraint {
                a,
                b,
                c,
                op: Op::MUL,
            }
        };
        res.push(new_constraint);
    }

    ConstraintSystem::new(&res)
}

#[cfg(test)]
mod test {
    use super::read_circom_constraint;
    use ark_bls12_381::Fr;

    #[test]
    fn test_read_circom_constraint() {
        let constraints =
            read_circom_constraint::<Fr>("src/circom/test_constraints.json".to_string());

        println!("Total constraints: {:#?}", constraints.total_constraint);
        println!("All constraints: {:#?}", constraints);
    }
}
