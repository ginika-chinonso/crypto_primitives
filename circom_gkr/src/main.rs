pub mod cli_functions;
pub mod utility_functions;
pub mod circom_adapter;

use clap::Parser;
use cli_functions::{Args, Commands};

pub fn main() {
    let args = Args::parse();

    let _ = match args.cmd {
        Commands::Prove {
            circom_file_path,
            public_input_path,
        } => Commands::prove(circom_file_path, public_input_path),
        Commands::Verify {
            circom_file_path,
            proof_path,
            public_input_path,
        } => Commands::verify_proof(circom_file_path, proof_path, public_input_path),
        Commands::ProveAndVerify {
            circom_file_path,
            public_input_path,
        } => Commands::prove_and_verify(circom_file_path, public_input_path),
        Commands::GenerateCircuit { circom_file_path } => {
            Commands::generate_circuit(circom_file_path)
        }
    };
}
