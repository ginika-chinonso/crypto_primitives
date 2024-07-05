#[cfg(test)]
pub mod test {
    use ark_bn254::Fq;
    use polynomials::multilinear_polynomial::eval_form::MLE;

    use crate::eval_form::{prover::Prover, verifier::Verifier};

    #[test]
    pub fn test_sumcheck_protocol() {
        let val = vec![
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(2),
            Fq::from(2),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(2),
            Fq::from(2),
        ];

        let v = vec![val; 4096 * 32].concat();

        dbg!(&v.len());

        let init_poly = MLE::new(&v);

        dbg!(&init_poly.num_of_vars);

        let proof = Prover::prove(init_poly);

        match Verifier::verify(proof) {
            Ok(val) => assert!(val, "Sumcheck failed"),
            Err(e) => panic!("Error occured: {}", e),
        }
    }
}
