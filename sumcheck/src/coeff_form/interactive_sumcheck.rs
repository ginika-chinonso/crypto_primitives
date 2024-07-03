// Sumcheck protocol

use ark_ff::PrimeField;
use polynomials::univariate_polynomial::UnivariatePolynomial;

pub struct Channel<F: PrimeField> {
    pub message: Vec<UnivariatePolynomial<F>>,
}

#[cfg(test)]
mod tests {
    use super::super::prover::Prover;
    use polynomials::{
        multilinear_polynomial::coef_form::{MultilinearMonomial, MultilinearPolynomial},
        univariate_polynomial::UnivariatePolynomial,
    };

    use ark_bls12_381::Fr;
    type Fq = Fr;

    // #[test]
    // fn test_sumcheck() {

    //     // At the begining of round, prover sends polynomial and the sum of evaluating the said poly at the boolean hypercube to the to verifier

    //     // Prover
    //     let poly = MultilinearPolynomial::<Fq>::interpolate(
    //         vec![
    //             Fq::from(3),
    //             Fq::from(5),
    //             Fq::from(2),
    //             Fq::from(9),
    //             Fq::from(11),
    //             Fq::from(16),
    //             Fq::from(1),
    //             Fq::from(5),
    //         ]
    //     );

    //     let mut num_of_vars = poly.terms[0].vars.len();
    //     let mut sum = Fq::from(0);
    //     for i in 0..2_usize.pow(num_of_vars as u32) {
    //         let eval_points = get_binary_string(i, num_of_vars);
    //         let eval_points_as_field : Vec<(usize, Fq)>= eval_points.chars().enumerate().map(|(a, b)| {
    //             if b == '0' {
    //                 (a, Fq::from(0))
    //             } else {
    //                 (a, Fq::from(1))
    //             }
    //         }).collect();
    //         let eval_res = &poly.evaluate(eval_points_as_field);
    //         sum += eval_res;
    //     };
    //     assert!(sum == Fq::from(52));

    //     // Prover generates a univariate polynomial by evaluating the polynomial at the boolean hypercube
    //     let mut round_0_poly = MultilinearPolynomial::<Fq>::new(vec![]);
    //     num_of_vars -= 1;
    //     let mut round = 0_usize;

    //     for i in 0..2_usize.pow(num_of_vars as u32) {
    //         let eval_points = get_binary_string(i, num_of_vars);
    //         let eval_points_as_field : Vec<(usize, Fq)>= eval_points.chars().enumerate().map(|(a, b)| {
    //             if b == '0' {
    //                 (a + round.clone() + 1, Fq::from(0))
    //             } else {
    //                 (a + round.clone() + 1, Fq::from(1))
    //             }
    //         }).collect();
    //         let eval_res = poly.clone().partial_eval(eval_points_as_field).relabel();
    //         round_0_poly = round_0_poly.add(eval_res.clone());
    //     };

    //     let round_0_univariate_poly = Polynomial::try_from(round_0_poly.relabel()).unwrap();

    //     // Prover sends univariate poly to verifier who then evaluates it at the boolean hypercube and checks that the sum equals the sum at the previous round
    //     let verifier_check = round_0_univariate_poly.evaluate(Fq::from(0)) + round_0_univariate_poly.evaluate(Fq::from(1));
    //     assert!(verifier_check == sum);

    //     // Round 1
    //     round += 1;

    //     // The verifier samples a random point and sends to the prover
    //     let round_1_challenge = Fq::from(5);
    //     let mut challenges = vec![];
    //     challenges.push(round_1_challenge);

    //     // Prover generates a univariate polynomial by evaluating the polynomial at the random point and the boolean hypercube
    //     let mut round_1_poly = MultilinearPolynomial::<Fq>::new(vec![]);
    //     num_of_vars -= 1;

    //     for i in 0..2_usize.pow(num_of_vars as u32) {
    //         let eval_points = get_binary_string(i, num_of_vars);
    //         let mut eval_points_as_field : Vec<(usize, Fq)>= eval_points.chars().enumerate().map(|(a, b)| {
    //             if b == '0' {
    //                 (a + round.clone() + 1, Fq::from(0))
    //             } else {
    //                 (a + round.clone() + 1, Fq::from(1))
    //             }
    //         }).collect();

    //         for i in 0..round {
    //             dbg!(&i);
    //             eval_points_as_field.push((i, challenges[i]));
    //         }
    //         let eval_res = poly.clone().partial_eval(eval_points_as_field).relabel();
    //         round_1_poly = round_1_poly.add(eval_res.clone());
    //     };

    //     let round_1_univariate_poly = Polynomial::try_from(round_1_poly.relabel()).unwrap();
    //     dbg!(&round_1_univariate_poly);

    //     // Prover sends univariate poly to verifier who then evaluates it at the boolean hypercube and checks that the sum equals the sum at the previous round
    //     let verifier_check = round_1_univariate_poly.evaluate(Fq::from(0)) + round_1_univariate_poly.evaluate(Fq::from(1));
    //     assert!(verifier_check == round_0_univariate_poly.evaluate(challenges[round - 1]));
    //     dbg!(&verifier_check);
    //     dbg!(round_0_univariate_poly.evaluate(challenges[round - 1]));

    // }

    #[test]
    fn test_sumcheck_prover() {
        let poly = MultilinearPolynomial::new(vec![
            MultilinearMonomial::new(Fq::from(2), vec![true, true, false]),
            MultilinearMonomial::new(Fq::from(3), vec![false, true, true]),
        ]);

        let prover = Prover::new(poly, Fq::from(10));

        let round_1_poly = prover.prove(&vec![]);
        assert!(round_1_poly == UnivariatePolynomial::new(vec![Fq::from(3), Fq::from(4)]));

        let round_2_poly = prover.prove(&vec![Fq::from(5)]);
        assert!(round_2_poly == UnivariatePolynomial::new(vec![Fq::from(0), Fq::from(23)]));

        let round_3_poly = prover.prove(&vec![Fq::from(5), Fq::from(9)]);
        assert!(round_3_poly == UnivariatePolynomial::new(vec![Fq::from(90), Fq::from(27)]));
    }

    // #[test]
    // fn test_sumcheck_verifier() {
    //     let poly = MultilinearPolynomial::new(vec![
    //         MultilinearMonomial::new(Fq::from(2), vec![true, true, false]),
    //         MultilinearMonomial::new(Fq::from(3), vec![false, true, true]),
    //     ]);

    //     let mut verifier = Verifier::new(poly, Fq::from(10));
    //     verifier.challenges.push(Fq::from(5));

    //     let round_1_verification = verifier.verify(Polynomial::new(vec![Fq::from(3), Fq::from(4)]));

    //     assert!(round_1_verification);

    //     verifier.challenges.push(Fq::from(9));
    //     let round_2_verification = verifier.verify(Polynomial::new(vec![Fq::from(0), Fq::from(23)]));

    //     assert!(round_2_verification);

    //     verifier.challenges.push(Fq::from(10));
    //     let round_3_verification = verifier.verify(Polynomial::new(vec![Fq::from(90), Fq::from(27)]));

    //     assert!(round_3_verification);
    // }
}
