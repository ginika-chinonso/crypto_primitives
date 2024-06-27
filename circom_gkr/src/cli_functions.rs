use std::error::Error;

use clap::{Parser, Subcommand};

use crate::utility_functions::{generate_circuit, prove, prove_and_verify, verify_proof};

#[derive(Debug, Parser)]
#[command(author, version, about, long_about = None)]
pub struct Args {
    #[command(subcommand)]
    pub cmd: Commands,
}

#[derive(Debug, Subcommand)]
pub enum Commands {
    Prove {
        #[clap(long)]
        circom_file_path: String,

        #[clap(long)]
        public_input_path: String,
    },

    Verify {
        #[clap(long)]
        circom_file_path: String,

        #[clap(long)]
        proof_path: String,

        #[clap(long)]
        public_input_path: String,
    },

    ProveAndVerify {
        #[clap(long)]
        circom_file_path: String,

        #[clap(long)]
        public_input_path: String,
    },

    GenerateCircuit {
        #[clap(long)]
        circom_file_path: String,
    },
}

impl Commands {
    pub fn prove(circom_file: String, witness_file: String) -> Result<(), Box<dyn Error>> {
        prove(circom_file, witness_file)?;
        Ok(())
    }

    pub fn generate_circuit(circom_file_path: String) -> Result<(), Box<dyn Error>> {
        generate_circuit(circom_file_path)?;
        Ok(())
    }

    pub fn verify_proof(
        circom_file_path: String,
        proof_path: String,
        public_input_path: String,
    ) -> Result<(), Box<dyn Error>> {
        verify_proof(circom_file_path, proof_path, public_input_path)?;
        Ok(())
    }

    pub fn prove_and_verify(
        circom_file_path: String,
        witness_file: String,
    ) -> Result<(), Box<dyn Error>> {
        prove_and_verify(circom_file_path, witness_file)?;
        Ok(())
    }
}
