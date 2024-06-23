use std::{error::Error, fs::File, io::{BufReader, BufWriter}, path::Path, process::{Command, Output}};
use ark_bn254::Fr;
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize, Read};
use gkr::{circuit::Circuit, protocol::{GKRProof, GKR}};
use serde_json::Value;
use crate::circom_adapter::{circom_read::read_circom_constraint, constraint_system::ConstraintSystem};


pub fn prove(circom_file_path: String, witness_file: String) -> Result<(), Box<dyn Error>> {

    let file_name: String = circom_file_path.split(".circom").collect();

    let _ = run_command(
        "circom",
        vec![
            circom_file_path.as_str(),
            "--O0",
            "--r1cs",
            "--sym",
            "--json",
            "--wasm",
        ],
        "Failed to compile circcom code"
    )?;

    let mut constraints: ConstraintSystem<Fr> = read_circom_constraint(format!("{}_constraints.json", file_name));

    constraints.simplify_constraint_system();

    let witness = generate_witness(file_name, witness_file, &constraints)?;
    
    let circuit: Circuit = constraints.into();

    println!("Starting to prove");

    let proof = GKR::prove(&circuit, &witness);

    println!("Proving finished");

    let file = File::create("proof.bin").unwrap();

    let writer = BufWriter::new(file);

    GKRProof::serialize_compressed(&proof, writer).unwrap();

    println!("Proof written to file: proof.bin");

    Ok(())
}

/// Takes the path to a circom program as input
/// Writes the circuit constraints generated to circuit.json
pub fn generate_circuit(circom_file_path: String) -> Result<(), Box<dyn Error>> {

    let file_name: String = circom_file_path.split(".circom").collect();

    let _ = run_command(
        "circom",
        vec![
            circom_file_path.as_str(),
            "--O0",
            "--r1cs",
            "--sym",
            "--json",
            "--wasm",
        ],
        "Failed to compile circcom code"
    )?;

    let mut constraints: ConstraintSystem<Fr> = read_circom_constraint(format!("{}_constraints.json", file_name));

    constraints.simplify_constraint_system();

    let circuit: Circuit = constraints.into();

    let file = File::create("circuit.json").expect("Failed to create circuit file");

    serde_json::to_writer(file, &circuit).unwrap();

    println!("Circuit generated and written to circuit.json file");

    Ok(())
}

pub fn verify_proof(
    circom_file_path: String,
    proof_path: String,
    public_input_path: String,
) -> Result<(), Box<dyn Error>> {

    let file_name: String = circom_file_path.split(".circom").collect();

    let _ = run_command(
        "circom",
        vec![
            circom_file_path.as_str(),
            "--O0",
            "--r1cs",
            "--sym",
            "--json",
            "--wasm",
        ],
        "Failed to compile circcom code"
    )?;

    let mut constraints: ConstraintSystem<Fr> = read_circom_constraint(format!("{}_constraints.json", file_name));

    constraints.simplify_constraint_system();

    let witness = generate_witness(file_name, public_input_path, &constraints)?;
    
    let circuit: Circuit = constraints.into();

    let proof_file = File::open(proof_path).expect("Failed to open proof file");

    let reader = BufReader::new(proof_file);

    let proof = GKRProof::deserialize_compressed(reader).expect("Unable to read proof");

    println!("Begining verification");

    let verify = GKR::verify(witness, proof.clone(), circuit).expect("Error during verification");

    assert!(verify, "Invalid GKR proof");

    // This assertion is necessary because 
    // the output of the circuit is always meant to be zero
    // This is not generalized but due to the implementation of the circuit
    assert!(proof.sumcheck_proofs[0].sum == Fr::from(0), "Invalid circuit eval");

    Ok(())
}

pub fn prove_and_verify(
    circom_file_path: String,
    witness_file: String,
) -> Result<(), Box<dyn Error>> {

    let file_name: String = circom_file_path.split(".circom").collect();

    let _ = run_command(
        "circom",
        vec![
            circom_file_path.as_str(),
            "--O0",
            "--r1cs",
            "--sym",
            "--json",
            "--wasm",
        ],
        "Failed to compile circcom code"
    )?;

    let mut constraints: ConstraintSystem<Fr> = read_circom_constraint(format!("{}_constraints.json", file_name));

    constraints.simplify_constraint_system();

    let witness = generate_witness(file_name, witness_file, &constraints)?;
    
    let circuit: Circuit = constraints.into();

    println!("Starting to prove");

    let proof = GKR::prove(&circuit, &witness);

    println!("Proving finished");

    let file = File::create("proof.bin").unwrap();

    let writer = BufWriter::new(file);

    GKRProof::serialize_compressed(&proof, writer).expect("Unable to write proof");

    let proof_file = File::open("proof.bin").expect("Failed to open file");

    let reader = BufReader::new(proof_file);

    let proof = GKRProof::deserialize_compressed(reader).expect("Unable to read proof");

    println!("Begining verification");

    let verify = GKR::verify(witness, proof.clone(), circuit)?;

    println!("Verification finished with: {verify}");

    assert!(verify, "Verification failed");

    // This assertion is necessary because 
    // the output of the circuit is always meant to be zero
    // This is not generalized but due to the implementation of the circuit
    assert!(proof.sumcheck_proofs[0].sum == Fr::from(0), "Invalid circuit eval");

    Ok(())
}

pub fn generate_witness(file_name: String, input_path: String, constraints: &ConstraintSystem<Fr>) -> Result<Vec<Fr>, Box<dyn Error>> {
    // node generate_witness.js multiplier2.wasm input.json witness.wtns

    let js_folder = format!("{}_js", file_name);

    let current_directory = Path::new(".");

    let new_path = current_directory.join(js_folder);

    let _ = run_command(
        "node",
        vec![
            format!("{}/generate_witness.js", new_path.display()).as_str(),
            format!("{}/{file_name}.wasm", new_path.display()).as_str(),
            format!("{}", input_path).as_str(),
            format!("{}/witness.wtns", current_directory.display()).as_str(),
        ],
        "Unable to run generate_witness.js"
    )?;

    let _ = run_command("snarkjs", vec!["wtns", "export", "json", "witness.wtns"], "Unable to run generate_witness.js")?;

    let _ = run_command("cargo", vec!["run", "generate-circuit", "--circom-file-path", "test.circom"], "Failed to generate circuit");

    // Read the witness
    let mut witness_file = File::open("witness.json").expect("Failed to open witness file");

    let mut witness_buffer = String::new();

    witness_file.read_to_string(&mut witness_buffer).unwrap();

    let witness_values: Vec<Value> =
        serde_json::from_str(&witness_buffer).expect("Unable to read public input");

    let mut witness = Vec::new();

    for val in witness_values {
        witness.push(Fr::from(
            serde_json::from_str::<i64>(val.as_str().unwrap()).unwrap(),
        ));
    }

    println!("Witness: {:?}", witness.len());
    println!("Total constraints: {:?}", constraints.total_constraint);
    println!("Total vars: {:?}", constraints.total_vars);
    println!("Constants: {:?}", constraints.constant_map.len());

    witness.extend(vec![Fr::from(0); constraints.constant_map.len()]);

    println!("{:?}", constraints.constant_map);

    println!("Witness: {:?}", witness);
    
    for (val, ind) in &constraints.constant_map {
        // if *ind == 0 {
        //     continue;
        // }
        println!("{}: {}", ind, val);
        witness[*ind] = Fr::from(*val);
    }
        
    
    println!("Witness: {:?}", witness);
        
    Ok(witness)

}

pub fn run_command(command: &str, arguments: Vec<&str>, error_message: &str) -> Result<Output, Box<dyn Error>> {
    match Command::new(&command).args(arguments).output() {
        Ok(res) => Ok(res),
        Err(_) => Err(error_message.into())
    }
}
