use ark_ff::{BigInteger, PrimeField};

use crate::{
    fiat_shamir_transcript::transcript::Transcript,
    gkr::utils::{eval_l, l_function, q_function},
    polynomials::{
        multilinear_poly::{MultilinearPolynomial, MultilinearPolynomialTrait},
        univariate_poly::UnivariatePolynomial,
    },
    sumcheck::non_interactive_sumcheck::{Sumcheck, SumcheckProof},
};

use super::{circuit::Circuit, eval_gate::EvalGate};

#[derive(Debug)]
pub struct GKRProof<F: PrimeField> {
    pub sumcheck_proofs: Vec<SumcheckProof<F, EvalGate<F>>>,
    pub output: MultilinearPolynomial<F>,
    pub q_funcs: Vec<UnivariatePolynomial<F>>,
}

impl<F: PrimeField> GKRProof<F> {
    pub fn new() -> Self {
        Self {
            sumcheck_proofs: vec![],
            output: MultilinearPolynomial::additive_identity(),
            q_funcs: vec![],
        }
    }
}

pub struct GKR {}

impl GKR {
    pub fn prove<F: PrimeField>(circuit: Circuit, input: Vec<F>) -> GKRProof<F> {
        let circuit_eval = circuit.evaluate(input);

        let mut transcript = Transcript::new();
        let mut gkr_proof = GKRProof::<F>::new();

        let output_gate_w_mle = circuit.w_mle(circuit_eval[0].clone());

        transcript.append(&output_gate_w_mle.to_bytes());

        let mut r = transcript.sample_n_field_elements(output_gate_w_mle.number_of_vars());

        let m_0_evalpoints = r.clone().into_iter().enumerate().collect();

        let mut m = output_gate_w_mle.clone().evaluate(&m_0_evalpoints);

        gkr_proof.output = output_gate_w_mle;

        for i in 1..circuit_eval.len() {
            dbg!(&i);
            dbg!(&circuit_eval.len());

            let [add_i_mle, mul_i_mle] = circuit.layer_mle(i - 1);

            let layer_w_mle = circuit.w_mle(circuit_eval[i].clone());

            let layer_eval_gate = EvalGate::new(
                &r,
                &vec![add_i_mle],
                &vec![mul_i_mle],
                &vec![layer_w_mle.clone()],
                &vec![layer_w_mle],
            );

            let round_sumcheck_proof = Sumcheck::prove(&layer_eval_gate, &m);
            // dbg!(&round_sumcheck_proof);
            // dbg!(&m);

            let (b, c) = round_sumcheck_proof
                .challenges
                .split_at(round_sumcheck_proof.challenges.len() / 2);

            let l = l_function(&b, &c);
            let q = q_function(&l, &circuit.w_mle(circuit_eval[i].clone())).unwrap();
            println!("Mid point");

            transcript.append(round_sumcheck_proof.to_bytes().as_slice());

            gkr_proof.sumcheck_proofs.push(round_sumcheck_proof);

            transcript.add_univariate_poly(&q);

            let r_star = transcript.sample_field_element();

            r = eval_l(&l, r_star);

            m = q.evaluate(r_star);

            gkr_proof.q_funcs.push(q);
        }

        println!("Done with proving");

        gkr_proof
    }

    pub fn verify<F: PrimeField>(
        input: Vec<F>,
        proof: GKRProof<F>,
        circuit: Circuit,
    ) -> Result<bool, String> {
        let mut transcript = Transcript::new();

        transcript.append(proof.output.to_bytes().as_slice());

        let mut r = transcript.sample_n_field_elements(proof.output.number_of_vars());

        let mut m = proof
            .output
            .evaluate(&r.clone().into_iter().enumerate().collect());

        // To check the provers claim, the verifier applies the sumcheck protocol to the polynomial
        // f(b, c) = add_i(a, b, c)(w_mle(b) + w_mle(c)) + add_i(a, b, c)(w_mle(b) + w_mle(c))

        for i in 0..proof.sumcheck_proofs.len() {
            transcript.append(proof.sumcheck_proofs[i].to_bytes().as_slice());
            transcript.add_univariate_poly(&proof.q_funcs[i]);

            dbg!(&m);
            dbg!(&proof.sumcheck_proofs[i].sum);
            if proof.sumcheck_proofs[i].sum != m {
                return Ok(false);
            }

            let (challenges, round_sum) =
                Sumcheck::verify_partial(&mut proof.sumcheck_proofs[i].clone())?;

            let challenges: Vec<F> = challenges
                .into_iter()
                .map(|(_, challenge)| challenge)
                .collect();

            let [add_i_mle, mul_i_mle] = circuit.layer_mle::<F>(i);

            let mut rbc = r.clone();

            rbc.extend(&challenges);

            let w_b = proof.q_funcs[i].evaluate(F::zero());

            let w_c = proof.q_funcs[i].evaluate(F::one());

            let (b, c) = challenges.split_at(challenges.len() / 2);

            // dbg!(&rbc);
            let rbc_eval_points: Vec<(usize, F)> = rbc.into_iter().enumerate().collect();

            let add_result = add_i_mle.evaluate(&rbc_eval_points) * (w_b + w_c);

            // dbg!(&add_i_mle.number_of_vars());
            // dbg!(&mul_i_mle.number_of_vars());
            // dbg!(&b);
            // dbg!(&c);
            // dbg!(&challenges);

            let mul_result = mul_i_mle.evaluate(&rbc_eval_points) * (w_b * w_c);

            let f_eval = add_result + mul_result;

            if f_eval != round_sum {
                return Ok(false);
            }

            let l_func = l_function(b, c);

            let r_star = transcript.sample_field_element();

            r = eval_l(&l_func, r_star);

            m = proof.q_funcs[i].evaluate(r_star);
        }

        let input_mle = MultilinearPolynomial::<F>::interpolate(input);

        let verifier_m = input_mle.evaluate(&r.into_iter().enumerate().collect());

        Ok(verifier_m == m)
    }
}

#[cfg(test)]
mod test {
    use ark_ff::{Fp64, MontBackend, MontConfig};

    use crate::{
        gkr::{
            circuit::{Circuit, Layer, Wire},
            eval_gate::EvalGate,
        },
        polynomials::multilinear_poly::{MultilinearPolynomial, MultilinearPolynomialTrait},
        sumcheck::non_interactive_sumcheck::Sumcheck,
    };

    use super::{GKRProof, GKR};

    #[derive(MontConfig)]
    #[modulus = "17"]
    #[generator = "3"]
    pub struct FqConfig;
    pub type Fq = Fp64<MontBackend<FqConfig, 1>>;

    fn create_circuit() -> Circuit {
        let layer0 = Layer::new(vec![], vec![Wire::new(0, 0, 1)]);
        let layer1 = Layer::new(vec![Wire::new(0, 0, 1)], vec![Wire::new(1, 1, 2)]);
        let layer2 = Layer::new(
            vec![Wire::new(0, 0, 1), Wire::new(2, 4, 5)],
            vec![Wire::new(1, 2, 3)],
        );

        let new_circuit = Circuit::new(vec![layer0, layer1, layer2]);

        new_circuit
    }

    #[test]
    fn test_gkr_proof() {
        let circuit = create_circuit();

        let circuit_input = vec![
            Fq::from(5),
            Fq::from(2),
            Fq::from(3),
            Fq::from(4),
            Fq::from(9),
            Fq::from(8),
        ];

        let gkr_proof: GKRProof<Fq> = GKR::prove(circuit.clone(), circuit_input.clone());

        // dbg!(&gkr_proof);

        let verifier = GKR::verify(circuit_input, gkr_proof, circuit).unwrap();

        assert!(verifier, "Proof not accepted");
    }

    #[test]
    fn test_sumcheck_eval() {
        let circuit = create_circuit();
        let circuit_input = vec![
            Fq::from(5),
            Fq::from(2),
            Fq::from(3),
            Fq::from(4),
            Fq::from(9),
            Fq::from(8),
        ];
        let circuit_eval = circuit.evaluate(circuit_input);
        let [add_i, mul_i] = circuit.layer_mle(1);
        let w_2 = circuit.w_mle(circuit_eval[2].clone());
        let gate_ext = EvalGate::new(
            &vec![Fq::from(0)],
            &vec![add_i],
            &vec![mul_i],
            &vec![w_2.clone()],
            &vec![w_2],
        );
        dbg!(gate_ext.number_of_vars());
        assert_eq!(
            Fq::from(2),
            gate_ext.evaluate(&vec![
                (0, Fq::from(0)),
                (1, Fq::from(0)),
                (2, Fq::from(0)),
                (3, Fq::from(1)),
            ])
        );

        let sumcheck_proof = Sumcheck::prove(&gate_ext, &Fq::from(2));
        assert!(Sumcheck::verify(sumcheck_proof, gate_ext).unwrap());
    }
}
