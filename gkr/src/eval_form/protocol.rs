use std::{error::Error, marker::PhantomData};

use ark_ff::PrimeField;
use fiat_shamir_transcript::Transcript;
use polynomials::multilinear_polynomial::traits::MultilinearPolynomialTrait;
use sumcheck::universal_mle::universal_mle::{Sumcheck, SumcheckProof};

use crate::eval_form::layer_poly::{evaluate_layer_mle_given_inputs, get_new_ops_mles};

use super::{
    circuit::{get_wmle, Circuit},
    layer_poly::LayerPoly,
};

#[derive(Debug)]
pub struct GKRProof<F: PrimeField> {
    pub output: Vec<F>,
    pub sumcheck_proofs: Vec<SumcheckProof<F>>,
    pub wbs: Vec<F>,
    pub wcs: Vec<F>,
}

impl<F: PrimeField> GKRProof<F> {
    pub fn new(
        output: Vec<F>,
        sumcheck_proofs: Vec<SumcheckProof<F>>,
        wbs: Vec<F>,
        wcs: Vec<F>,
    ) -> Self {
        Self {
            output,
            sumcheck_proofs,
            wbs,
            wcs,
        }
    }
}

pub struct GKR<F> {
    _marker: PhantomData<F>,
}

impl<F: PrimeField> GKR<F> {
    pub fn prove(circuit: &Circuit<F>, inputs: &Vec<F>) -> GKRProof<F> {
        let circuit_evaluation = circuit.evaluate(inputs);

        let mut rounds_sumcheck_proof = vec![];

        let mut wbs = vec![];

        let mut wcs = vec![];

        let mut transcript = Transcript::new();

        let r: Vec<F> =
            transcript.sample_n_field_elements(circuit.layers[0].layer_output_num_of_vars());

        let claimed_sum = &circuit_evaluation[0];

        let layer_ops_mle = circuit.get_layer_ops_mles(0);

        let layer_w_mle = get_wmle(&circuit_evaluation[1]);

        let layer_poly = LayerPoly::new(&r, &layer_ops_mle, &layer_w_mle);

        let layer_sumcheck_proof = Sumcheck::prove(layer_poly.poly);

        rounds_sumcheck_proof.push(layer_sumcheck_proof.0);

        let (r_b0, r_c0) = layer_sumcheck_proof
            .1
            .split_at(layer_sumcheck_proof.1.len() / 2);

        let mut r_b = r_b0.to_vec();
        let mut r_c = r_c0.to_vec();

        wbs.push(
            layer_w_mle.evaluate(
                &r_b.iter()
                    .enumerate()
                    .map(|(ind, val)| (ind + 1, *val))
                    .collect(),
            ),
        );

        wcs.push(
            layer_w_mle.evaluate(
                &r_c.iter()
                    .enumerate()
                    .map(|(ind, val)| (ind + 1, *val))
                    .collect(),
            ),
        );

        for i in 1..circuit.depth {
            let layer_i_ops_mle = circuit.get_layer_ops_mles(i);

            let alpha_n_beta = transcript.sample_n_field_elements(2);

            let new_ops_mles = get_new_ops_mles(&layer_i_ops_mle, &alpha_n_beta, &r_b, &r_c);

            let layer_i_plus_1_w_mle = get_wmle(&circuit_evaluation[i + 1]);

            let new_layer_poly = LayerPoly::new(&vec![], &new_ops_mles, &layer_i_plus_1_w_mle);

            let layer_sumcheck_proof = Sumcheck::prove(new_layer_poly.poly);

            rounds_sumcheck_proof.push(layer_sumcheck_proof.0);

            let (r_b0, r_c0) = layer_sumcheck_proof
                .1
                .split_at(layer_sumcheck_proof.1.len() / 2);

            r_b = r_b0.to_vec();

            r_c = r_c0.to_vec();

            wbs.push(
                layer_i_plus_1_w_mle.evaluate(
                    &r_b.iter()
                        .enumerate()
                        .map(|(ind, val)| (ind + 1, *val))
                        .collect(),
                ),
            );

            wcs.push(
                layer_i_plus_1_w_mle.evaluate(
                    &r_c.iter()
                        .enumerate()
                        .map(|(ind, val)| (ind + 1, *val))
                        .collect(),
                ),
            );
        }

        GKRProof::new(claimed_sum.clone(), rounds_sumcheck_proof, wbs, wcs)
    }

    pub fn verify(
        circuit: &Circuit<F>,
        inputs: &Vec<F>,
        proof: &GKRProof<F>,
    ) -> Result<bool, Box<dyn Error>> {
        let mut transcript = Transcript::new();

        let mut r: Vec<F> =
            transcript.sample_n_field_elements(circuit.layers[0].layer_output_num_of_vars());

        let mut claimed_sum = get_wmle(&proof.output).evaluate(
            &r.iter()
                .enumerate()
                .map(|(ind, val)| (ind + 1, *val))
                .collect(),
        );

        assert_eq!(
            claimed_sum, proof.sumcheck_proofs[0].claimed_sum,
            "Claimed sum does not match for 0"
        );

        let mut partial_verification = Sumcheck::verify_partial(&proof.sumcheck_proofs[0]);

        let (mut challenges, mut round_sum): (Vec<F>, F) = match partial_verification {
            Ok(val) => (
                val.0.iter().map(|(_, challenge)| *challenge).collect(),
                val.2,
            ),
            Err(err) => panic!("Partial verification failed: {}", err),
        };

        r.extend(challenges.clone());

        let (r_b0, r_c0) = challenges.split_at(challenges.len() / 2);

        let mut r_b = r_b0.to_vec();

        let mut r_c = r_c0.to_vec();

        let f_b_c_eval = evaluate_layer_mle_given_inputs(
            &circuit.get_layer_ops_mles(0),
            &r,
            &proof.wbs[0],
            &proof.wcs[0],
        );

        assert_eq!(f_b_c_eval, round_sum, "Round sum does not match");

        let mut alpha_n_beta: Vec<F> = transcript.sample_n_field_elements(2);

        claimed_sum = alpha_n_beta[0] * &proof.wbs[0] + alpha_n_beta[1] * proof.wcs[0];

        let mut new_layer_ops_mle = vec![];

        for i in 1..(proof.sumcheck_proofs.len() - 1) {
            assert_eq!(
                claimed_sum, proof.sumcheck_proofs[i].claimed_sum,
                "Claimed sum does not match for {i}"
            );

            new_layer_ops_mle =
                get_new_ops_mles(&circuit.get_layer_ops_mles(i), &alpha_n_beta, &r_b, &r_c);

            partial_verification = Sumcheck::verify_partial(&proof.sumcheck_proofs[i]);

            (challenges, round_sum) = match partial_verification {
                Ok(val) => (
                    val.0.iter().map(|(_, challenge)| *challenge).collect(),
                    val.2,
                ),
                Err(err) => panic!("Partial verification failed: {}", err),
            };

            let (r_b0, r_c0) = challenges.split_at(challenges.len() / 2);

            r_b = r_b0.to_vec();

            r_c = r_c0.to_vec();

            let f_b_c_eval = evaluate_layer_mle_given_inputs(
                &new_layer_ops_mle,
                &challenges,
                &proof.wbs[i],
                &proof.wcs[i],
            );

            assert_eq!(f_b_c_eval, round_sum, "Round sum does not match");

            alpha_n_beta = transcript.sample_n_field_elements(2);

            claimed_sum = alpha_n_beta[0] * &proof.wbs[i] + alpha_n_beta[1] * proof.wcs[i];
        }

        new_layer_ops_mle = get_new_ops_mles(
            &circuit.get_layer_ops_mles(circuit.depth - 1),
            &alpha_n_beta,
            &r_b,
            &r_c,
        );

        let partial_verification =
            Sumcheck::verify_partial(&proof.sumcheck_proofs[proof.sumcheck_proofs.len() - 1]);

        (challenges, round_sum) = match partial_verification {
            Ok(val) => (
                val.0.iter().map(|(_, challenge)| *challenge).collect(),
                val.2,
            ),
            Err(err) => panic!("Partial verification failed: {}", err),
        };

        claimed_sum = evaluate_layer_mle_given_inputs(
            &new_layer_ops_mle,
            &challenges,
            &proof.wbs[proof.wbs.len() - 1],
            &proof.wcs[proof.wcs.len() - 1],
        );

        assert_eq!(
            round_sum, claimed_sum,
            "Round sum and claimed sum does not match"
        );

        let input_wle = get_wmle(&inputs);

        let layer_poly = LayerPoly::new(&vec![], &new_layer_ops_mle, &input_wle);

        let input_layer_poly_eval = layer_poly.evaluate(
            &challenges
                .iter()
                .enumerate()
                .map(|(ind, val)| (ind + 1, *val))
                .collect(),
        );

        assert_eq!(claimed_sum, input_layer_poly_eval, "Invalid claimed sum");

        Ok(true)
    }
}

#[cfg(test)]
pub mod tests {

    use crate::eval_form::circuit::{Circuit, Gate, Layer};
    use ark_bn254::Fq;
    use ark_ff::PrimeField;
    use polynomials::multilinear_polynomial::universal_mle::ops::Ops;

    use super::GKR;

    pub fn create_circuit<F: PrimeField>() -> Circuit<F> {
        let ops = Ops::new();

        let layer_4 = Layer::new(vec![
            Gate::new(0, 0, 1, &ops.add),
            Gate::new(1, 2, 3, &ops.mul),
            Gate::new(2, 4, 5, &ops.add),
            Gate::new(3, 6, 7, &ops.mul),
            Gate::new(4, 8, 9, &ops.mul),
            Gate::new(5, 10, 11, &ops.add),
            Gate::new(6, 12, 13, &ops.mul),
            Gate::new(7, 14, 15, &ops.add),
            Gate::new(8, 16, 17, &ops.mul),
            Gate::new(9, 18, 19, &ops.mul),
            Gate::new(10, 20, 21, &ops.add),
            Gate::new(11, 22, 23, &ops.mul),
        ]);

        let layer_3 = Layer::new(vec![
            Gate::new(0, 0, 1, &ops.add),
            Gate::new(1, 2, 3, &ops.mul),
            Gate::new(2, 4, 5, &ops.add),
            Gate::new(3, 6, 7, &ops.mul),
            Gate::new(4, 8, 9, &ops.mul),
            Gate::new(5, 10, 11, &ops.add),
        ]);

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

        Circuit::new(vec![layer_0, layer_1, layer_2, layer_3, layer_4])
    }

    #[test]
    pub fn test_proving() {
        let circuit = create_circuit();

        let inputs = vec![
            Fq::from(1),
            Fq::from(2),
            Fq::from(3),
            Fq::from(4),
            Fq::from(5),
            Fq::from(6),
            Fq::from(7),
            Fq::from(8),
            Fq::from(9),
            Fq::from(10),
            Fq::from(11),
            Fq::from(12),
            Fq::from(13),
            Fq::from(14),
            Fq::from(15),
            Fq::from(16),
            Fq::from(17),
            Fq::from(18),
            Fq::from(19),
            Fq::from(20),
            Fq::from(21),
            Fq::from(22),
            Fq::from(23),
            Fq::from(24),
        ];

        let proof = GKR::prove(&circuit, &inputs);

        let verify = GKR::verify(&circuit, &inputs, &proof);

        assert!(verify.unwrap(), "Invalid proof");
    }
}
