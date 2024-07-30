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

        let init_poly = MLE::new(&v);

        let proof = Prover::prove_sumcheck(init_poly);

        match Verifier::verify_sumcheck(proof) {
            Ok(val) => assert!(val, "Sumcheck failed"),
            Err(e) => panic!("Error occured: {}", e),
        }
    }

    #[test]
    pub fn test_sum_of_product_protocol() {
        let f_poly = MLE::new(&vec![
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
        ]);

        let g_poly = MLE::new(&vec![
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(3),
            Fq::from(3),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(3),
            Fq::from(3),
        ]);

        let proof = Prover::prove_sum_of_product(&f_poly, &g_poly);

        let verify = Verifier::verify_sumcheck_product(proof);

        assert!(verify.unwrap(), "Invalid sumcheck proof");
    }
}
