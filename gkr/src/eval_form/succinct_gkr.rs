use std::{error::Error, marker::PhantomData};

use ark_ec::pairing::Pairing;
use ark_ff::{field_hashers::DefaultFieldHasher, BigInteger, PrimeField};
use fiat_shamir_transcript::Transcript;
use multilinear_kzg::{MultilinearPCS, MultilinearPCSProof};
use polynomials::multilinear_polynomial::{eval_form::MLE, traits::MultilinearPolynomialTrait};
use sha2::Sha256;
use sumcheck::universal_mle::universal_mle::{Sumcheck, SumcheckProof};

use crate::eval_form::layer_poly::{evaluate_layer_mle_given_inputs, get_new_ops_mles};

use super::{
    circuit::{get_wmle, Circuit},
    layer_poly::LayerPoly,
};

#[derive(Debug)]
pub struct GKRProof<P: Pairing> {
    pub output: Vec<P::ScalarField>,
    pub sumcheck_proofs: Vec<SumcheckProof<P::ScalarField>>,
    pub wbs: Vec<P::ScalarField>,
    pub wcs: Vec<P::ScalarField>,
    pub input_poly_commitment: P::G1,
    pub input_rb_opening: MultilinearPCSProof<P>,
    pub input_rc_opening: MultilinearPCSProof<P>,
}

impl<P: Pairing> GKRProof<P> {
    pub fn new(
        output: Vec<P::ScalarField>,
        sumcheck_proofs: Vec<SumcheckProof<P::ScalarField>>,
        wbs: Vec<P::ScalarField>,
        wcs: Vec<P::ScalarField>,
        input_poly_commitment: P::G1,
        input_rb_opening: MultilinearPCSProof<P>,
        input_rc_opening: MultilinearPCSProof<P>,
    ) -> Self {
        Self {
            output,
            sumcheck_proofs,
            wbs,
            wcs,
            input_poly_commitment,
            input_rb_opening,
            input_rc_opening,
        }
    }
}

pub struct GKR<P> {
    _pairing: PhantomData<P>,
}

impl<P: Pairing> GKR<P> {
    pub fn prove(
        circuit: &Circuit<P::ScalarField>,
        inputs: &Vec<P::ScalarField>,
        multilinear_pcs: &MultilinearPCS<P, DefaultFieldHasher<Sha256>>,
    ) -> GKRProof<P> {
        let circuit_evaluation = circuit.evaluate(inputs);

        let input_mle = MLE::new(&inputs);

        let input_poly_commitment = multilinear_pcs.commit(&input_mle);

        let mut rounds_sumcheck_proof = vec![];

        let mut wbs = vec![];

        let mut wcs = vec![];

        let mut transcript = Transcript::new();

        transcript.append(&circuit_evaluation[0][0].into_bigint().to_bytes_be());
        transcript.append(&input_poly_commitment.to_string().into_bytes());

        let r: Vec<P::ScalarField> =
            transcript.sample_n_field_elements(circuit.layers[0].layer_output_num_of_vars());

        let claimed_sum = &circuit_evaluation[0];

        let layer_ops_mle = circuit.get_layer_ops_mles(0);

        let layer_w_mle = get_wmle(&circuit_evaluation[1]);

        let layer_poly = LayerPoly::new(&r, &layer_ops_mle, &layer_w_mle);

        let layer_sumcheck_proof = Sumcheck::prove(layer_poly.poly);

        transcript.append(&layer_sumcheck_proof.0.to_bytes());

        rounds_sumcheck_proof.push(layer_sumcheck_proof.0);

        let (r_b0, r_c0) = layer_sumcheck_proof
            .1
            .split_at(layer_sumcheck_proof.1.len() / 2);

        let mut r_b = r_b0.to_vec();
        let mut r_c = r_c0.to_vec();

        let mut wb = layer_w_mle.evaluate(
            &r_b.iter()
                .enumerate()
                .map(|(ind, val)| (ind + 1, *val))
                .collect(),
        );

        let mut wc = layer_w_mle.evaluate(
            &r_c.iter()
                .enumerate()
                .map(|(ind, val)| (ind + 1, *val))
                .collect(),
        );

        transcript.append(&wb.into_bigint().to_bytes_be());
        transcript.append(&wc.into_bigint().to_bytes_be());

        wbs.push(wb);

        wcs.push(wc);

        for i in 1..circuit.depth {
            let layer_i_ops_mle = circuit.get_layer_ops_mles(i);

            let alpha_n_beta = transcript.sample_n_field_elements(2);

            let new_ops_mles = get_new_ops_mles(&layer_i_ops_mle, &alpha_n_beta, &r_b, &r_c);

            let layer_i_plus_1_w_mle = get_wmle(&circuit_evaluation[i + 1]);

            let new_layer_poly = LayerPoly::new(&vec![], &new_ops_mles, &layer_i_plus_1_w_mle);

            let layer_sumcheck_proof = Sumcheck::prove(new_layer_poly.poly);

            transcript.append(&layer_sumcheck_proof.0.to_bytes());

            rounds_sumcheck_proof.push(layer_sumcheck_proof.0);

            let (r_b0, r_c0) = layer_sumcheck_proof
                .1
                .split_at(layer_sumcheck_proof.1.len() / 2);

            r_b = r_b0.to_vec();

            r_c = r_c0.to_vec();

            wb = layer_i_plus_1_w_mle.evaluate(
                &r_b.iter()
                    .enumerate()
                    .map(|(ind, val)| (ind + 1, *val))
                    .collect(),
            );

            wc = layer_i_plus_1_w_mle.evaluate(
                &r_c.iter()
                    .enumerate()
                    .map(|(ind, val)| (ind + 1, *val))
                    .collect(),
            );

            transcript.append(&wb.into_bigint().to_bytes_be());
            transcript.append(&wc.into_bigint().to_bytes_be());

            wbs.push(wb);

            wcs.push(wc);
        }

        let input_rb_opening = multilinear_pcs.open(
            &input_mle,
            &r_b.iter()
                .enumerate()
                .map(|(ind, val)| (ind + 1, *val))
                .collect(),
        );

        let input_rc_opening = multilinear_pcs.open(
            &input_mle,
            &r_c.iter()
                .enumerate()
                .map(|(ind, val)| (ind + 1, *val))
                .collect(),
        );

        GKRProof::new(
            claimed_sum.clone(),
            rounds_sumcheck_proof,
            wbs,
            wcs,
            input_poly_commitment,
            input_rb_opening,
            input_rc_opening,
        )
    }

    pub fn verify(
        circuit: &Circuit<P::ScalarField>,
        proof: &GKRProof<P>,
        multilinear_pcs: &MultilinearPCS<P, DefaultFieldHasher<Sha256>>,
    ) -> Result<bool, Box<dyn Error>> {
        let mut transcript = Transcript::new();

        transcript.append(&proof.output[0].into_bigint().to_bytes_be());
        transcript.append(&proof.input_poly_commitment.to_string().into_bytes());

        let mut r: Vec<P::ScalarField> =
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

        let (mut challenges, mut round_sum): (Vec<P::ScalarField>, P::ScalarField) =
            match partial_verification {
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

        transcript.append(&proof.sumcheck_proofs[0].to_bytes());
        transcript.append(&proof.wbs[0].into_bigint().to_bytes_be());
        transcript.append(&proof.wcs[0].into_bigint().to_bytes_be());

        let mut alpha_n_beta: Vec<P::ScalarField> = transcript.sample_n_field_elements(2);

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

            transcript.append(&proof.sumcheck_proofs[i].to_bytes());
            transcript.append(&proof.wbs[i].into_bigint().to_bytes_be());
            transcript.append(&proof.wcs[i].into_bigint().to_bytes_be());

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

        let (r_b0, r_c0) = challenges.split_at(challenges.len() / 2);

        r_b = r_b0.to_vec();

        r_c = r_c0.to_vec();

        let wb = if multilinear_pcs.verify(
            proof.input_poly_commitment,
            r_b.to_vec(),
            proof.input_rb_opening.clone(),
        ) {
            proof.input_rb_opening.evaluation
        } else {
            panic!("Invalid wb opening")
        };

        let wc = if multilinear_pcs.verify(
            proof.input_poly_commitment,
            r_c.to_vec(),
            proof.input_rc_opening.clone(),
        ) {
            proof.input_rc_opening.evaluation
        } else {
            panic!("Invalid wc opening");
        };

        let input_layer_poly_eval =
            evaluate_layer_mle_given_inputs(&new_layer_ops_mle, &challenges, &wb, &wc);

        assert_eq!(claimed_sum, input_layer_poly_eval, "Invalid claimed sum");

        Ok(true)
    }
}

#[cfg(test)]
pub mod tests {

    use crate::eval_form::circuit::{Circuit, Gate, Layer};
    use ark_bls12_381::{Bls12_381, Fr};
    use ark_ff::{field_hashers::DefaultFieldHasher, PrimeField};
    use multilinear_kzg::MultilinearPCS;
    use polynomials::multilinear_polynomial::universal_mle::ops::Ops;
    use sha2::Sha256;

    use super::GKR;

    type Fq = Fr;

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
    pub fn test_succinct_gkr() {
        let circuit = create_circuit();

        let mut inputs = vec![
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

        inputs.extend(vec![
            Fq::from(0);
            inputs.len().next_power_of_two() - inputs.len()
        ]);

        let taus = vec![
            "contrib for a".as_bytes().to_vec(),
            "contrib for b".as_bytes().to_vec(),
            "contrib for c".as_bytes().to_vec(),
            "contrib for d".as_bytes().to_vec(),
            "contrib for e".as_bytes().to_vec(),
        ];

        let multilinear_pcs: MultilinearPCS<Bls12_381, DefaultFieldHasher<Sha256>> =
            MultilinearPCS::instantiate("initialize multilinear pcs".as_bytes(), taus);

        let proof = GKR::prove(&circuit, &inputs, &multilinear_pcs);

        let verify = GKR::verify(&circuit, &proof, &multilinear_pcs);

        assert!(verify.unwrap(), "Invalid proof");
    }
}
