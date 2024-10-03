use ark_ff::PrimeField;
use polynomials::{
    multilinear_polynomial::{
        eval_form::MLE,
        traits::MultilinearPolynomialTrait,
        universal_mle::{ops::Ops, universal_mle::UniversalMLE},
    },
    univariate_polynomial::UnivariatePolynomial,
};

// f(a,rb,rc) = addi(a,rb,rc)(w(rb) + w(rc)) + multi(a,rb,rc)(w(rb)*w(rc))
// we need to reduce the verification of the two claims w(rb) and w(rc) to the verification of one claim
// we can do that by using a random linear combination of w(rb) and w(rc)
// we can as well prove that the w(rb) and w(rc) can be gotten from the layer beneath it
// by proving that the random linear combination can be gotten from the layer beneath
// where B = w(b) and C = w(c)
// alpha * B + beta * c = alpha * (addi+1(rb,b,c)(w(b) + w(c)) + multi(rb,b,c)(w(b)*w(c))) + beta * (addi+1(rc,b,c)(w(b)+w(c)) + multi(rc,b,c)(w(b)*w(c)))
// = (alpha * addi+1(rb,b,c)(w(b) + w(c)) + alpha * multi(rb,b,c)(w(b)*w(c))) + (beta * addi+1(rc,b,c)(w(b)+w(c)) + beta * multi(rc,b,c)(w(b)*w(c)))
// = (alpha * addi+1(rb,b,c)(w(b) + w(c)) + (beta * addi+1(rc,b,c)(w(b)+w(c)) + (alpha * multi(rb,b,c)(w(b)*w(c)) + (beta * multi(rc,b,c)(w(b)*w(c))
// = ((alpha * addi+1(rb,b,c)) + (beta * addi+1(rc,b,c)) (w(b)+w(c))) + (((alpha * multi(rb,b,c) + (beta * multi(rc,b,c)) (w(b)*w(c)))

#[derive(Debug, Clone)]
pub struct LayerPoly<F: PrimeField> {
    pub poly: UniversalMLE<F>,
}

impl<F: PrimeField> LayerPoly<F> {
    pub fn new(
        r: &Vec<F>,
        ops_mle: &Vec<UniversalMLE<F>>,
        w_mle: &UniversalMLE<F>,
    ) -> LayerPoly<F> {
        let r: Vec<(usize, F)> = r
            .into_iter()
            .enumerate()
            .map(|(index, val)| (index + 1, *val))
            .collect();

        let ops_mle_mul_wb_op_wc = ops_mle
            .iter()
            .map(|op_mle| op_mle.partial_eval(&r))
            .map(|op_mle| {
                let op = match &op_mle {
                    UniversalMLE::MLES(_, op) => op.clone(),
                    UniversalMLE::Value(_) => {
                        panic!("Gate doesnt have an operation")
                    }
                };

                let mut indices_for_wb = ((w_mle.number_of_vars() + 1)
                    ..=(2 * w_mle.number_of_vars()))
                    .collect::<Vec<usize>>();

                let mut indices_for_wc = (1..=(w_mle.number_of_vars())).collect();

                let wb_op_wc = UniversalMLE::new(
                    vec![
                        w_mle.clone().add_variable_at_index(&mut indices_for_wb),
                        w_mle.clone().add_variable_at_index(&mut indices_for_wc),
                    ],
                    op.clone(),
                );

                UniversalMLE::MLES(vec![op_mle, wb_op_wc], Ops::new().mul)
            })
            .collect();

        let poly = UniversalMLE::new(ops_mle_mul_wb_op_wc, Ops::new().add);

        Self { poly }
    }

    pub fn skip_one_and_sum_over_the_boolean_hypercube(&self) -> UnivariatePolynomial<F> {
        self.poly.skip_one_and_sum_over_the_boolean_hypercube()
    }
}

impl<F: PrimeField> MultilinearPolynomialTrait<F> for LayerPoly<F> {
    fn partial_eval(&self, x: &Vec<(usize, F)>) -> Self {
        LayerPoly {
            poly: self.poly.partial_eval(x),
        }
    }

    fn evaluate(&self, x: &Vec<(usize, F)>) -> F {
        self.poly.evaluate(x)
    }

    fn number_of_vars(&self) -> usize {
        self.poly.number_of_vars()
    }

    fn to_bytes(&self) -> Vec<u8> {
        self.poly.to_bytes()
    }

    fn relabel(&self) -> Self {
        todo!()
    }

    fn additive_identity() -> Self {
        LayerPoly::new(&vec![], &vec![], &UniversalMLE::additive_identity())
    }

    fn sum_over_the_boolean_hypercube(&self) -> F {
        self.poly.sum_over_the_boolean_hypercube()
    }

    fn to_univariate(
        &self,
    ) -> Result<polynomials::univariate_polynomial::UnivariatePolynomial<F>, String> {
        todo!()
    }
}

pub fn get_new_ops_mles<F: PrimeField>(
    ops_mle: &Vec<UniversalMLE<F>>,
    alpha_n_beta: &Vec<F>,
    r_b: &Vec<F>,
    r_c: &Vec<F>,
) -> Vec<UniversalMLE<F>> {
    ops_mle
        .iter()
        .map(|mle| {
            let mle_op = match mle {
                UniversalMLE::MLES(_, op) => op,
                UniversalMLE::Value(_) => panic!("MLE has no Operation"),
            };
            let alpha_mul_op_mle = mle
                .clone()
                .scalar_mul(alpha_n_beta[0])
                .partial_eval(
                    &r_b.iter()
                        .enumerate()
                        .map(|(ind, val)| (ind + 1, *val))
                        .collect(),
                )
                .to_eval();
            let beta_mul_op_mle = mle
                .clone()
                .scalar_mul(alpha_n_beta[1])
                .partial_eval(
                    &r_c.iter()
                        .enumerate()
                        .map(|(ind, val)| (ind + 1, *val))
                        .collect(),
                )
                .to_eval();
            let new_values = alpha_mul_op_mle
                .iter()
                .zip(beta_mul_op_mle)
                .map(|(lhs, rhs)| *lhs + rhs)
                .collect();

            UniversalMLE::MLES(
                vec![UniversalMLE::Value(MLE::new(&new_values))],
                mle_op.clone(),
            )
        })
        .collect()
}

pub fn evaluate_layer_mle_given_inputs<F: PrimeField>(
    layer_ops_mle: &Vec<UniversalMLE<F>>,
    challenges: &Vec<F>,
    wb: &F,
    wc: &F,
) -> F {
    layer_ops_mle.iter().fold(F::zero(), |mut init, mle| {
        let val = match mle {
            UniversalMLE::MLES(mles, op) => {
                let val1 = mles.iter().fold(F::zero(), |mut init, mle| {
                    init += mle.evaluate(
                        &challenges
                            .iter()
                            .enumerate()
                            .map(|(ind, val)| (ind + 1, *val))
                            .collect(),
                    ) * (op.op)(*wb, *wc);
                    init
                });
                val1
            }
            UniversalMLE::Value(_) => F::zero(),
        };
        init += val;
        init
    })
}

#[cfg(test)]
pub mod tests {
    use ark_bn254::Fq;
    use ark_ff::PrimeField;
    use polynomials::multilinear_polynomial::{
        traits::MultilinearPolynomialTrait,
        universal_mle::{ops::Ops, universal_mle::UniversalMLE},
    };

    use crate::eval_form::{
        circuit::{get_wmle, Circuit, Gate, Layer},
        layer_poly::get_new_ops_mles,
    };

    use super::{evaluate_layer_mle_given_inputs, LayerPoly};

    pub fn create_circuit<F: PrimeField>() -> Circuit<F> {
        let ops = Ops::new();

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

        Circuit::new(vec![layer_0, layer_1, layer_2])
    }

    pub fn create_circuit_1<F: PrimeField>() -> Circuit<F> {
        let ops = Ops::new();

        let layer_2 = Layer::new(vec![
            Gate::new(0, 0, 1, &ops.add),
            Gate::new(1, 2, 3, &ops.mul),
            Gate::new(2, 4, 5, &ops.add),
        ]);

        let layer_1 = Layer::new(vec![
            Gate::new(0, 0, 1, &ops.mul),
            Gate::new(1, 1, 2, &ops.mul),
        ]);

        let layer_0 = Layer::new(vec![Gate::new(0, 0, 1, &ops.mul)]);

        Circuit::new(vec![layer_0, layer_1, layer_2])
    }

    #[test]
    pub fn test_create_layer_op_poly() {
        let circuit = create_circuit();

        let inputs = vec![
            Fq::from(1),
            Fq::from(2),
            Fq::from(3),
            Fq::from(4),
            Fq::from(5),
            Fq::from(6),
        ];

        let circuit_evaluation = circuit.evaluate(&inputs);

        // Create op_add poly for layer 2
        let layer_op_mle = circuit.get_layer_op_mle(2, Ops::new().add);

        assert_eq!(
            layer_op_mle.sum_over_the_boolean_hypercube(),
            Fq::from(2),
            "Wrong number of add gates"
        );

        // Layer w_mle for Layer 2 evaluations
        let layer_w_mle = get_wmle(&circuit_evaluation[3]);

        // Layer op_add mle for output gate 00
        let layer_2_poly = LayerPoly::new(
            &vec![Fq::from(0), Fq::from(0)],
            &vec![layer_op_mle.clone()],
            &layer_w_mle,
        );
        assert_eq!(
            layer_2_poly.evaluate(&vec![
                (1, Fq::from(0)),
                (2, Fq::from(0)),
                (3, Fq::from(0)),
                (4, Fq::from(0)),
                (5, Fq::from(0)),
                (6, Fq::from(1))
            ]),
            Fq::from(3),
            "Wrong circuit evaluation"
        );

        // Layer op_add mle for output gate 10
        let layer_2_poly = LayerPoly::new(
            &vec![Fq::from(1), Fq::from(0)],
            &vec![layer_op_mle.clone()],
            &layer_w_mle,
        );
        assert_eq!(
            layer_2_poly.evaluate(&vec![
                (1, Fq::from(1)),
                (2, Fq::from(0)),
                (3, Fq::from(0)),
                (4, Fq::from(1)),
                (5, Fq::from(0)),
                (6, Fq::from(1))
            ]),
            Fq::from(11),
            "Wrong circuit evaluation"
        );

        // Create op_mul poly for layer 2
        let layer_op_mle = circuit.get_layer_op_mle(2, Ops::new().mul);

        // Layer op_mul mle for output gate 01
        let layer_2_poly = LayerPoly::new(
            &vec![Fq::from(0), Fq::from(1)],
            &vec![layer_op_mle],
            &layer_w_mle,
        );
        assert_eq!(
            layer_2_poly.evaluate(&vec![
                (1, Fq::from(0)),
                (2, Fq::from(1)),
                (3, Fq::from(0)),
                (4, Fq::from(0)),
                (5, Fq::from(1)),
                (6, Fq::from(1))
            ]),
            Fq::from(12),
            "Wrong circuit evaluation"
        );
    }

    #[test]
    pub fn create_layer_poly() {
        let circuit = create_circuit();

        let inputs = vec![
            Fq::from(1),
            Fq::from(2),
            Fq::from(3),
            Fq::from(4),
            Fq::from(5),
            Fq::from(6),
        ];

        let circuit_evaluation = circuit.evaluate(&inputs);

        // Create layer 2 poly
        let layer_op_mle = circuit.get_layer_ops_mles(2);

        // Layer w_mle for Layer 2 evaluations
        let layer_w_mle = get_wmle(&circuit_evaluation[3]);

        // Layer poly checking output gate 00
        assert_eq!(
            LayerPoly::new(&vec![Fq::from(0), Fq::from(0)], &layer_op_mle, &layer_w_mle,).evaluate(
                &vec![
                    (1, Fq::from(0)),
                    (2, Fq::from(0)),
                    (3, Fq::from(0)),
                    (4, Fq::from(0)),
                    (5, Fq::from(0)),
                    (6, Fq::from(1))
                ]
            ),
            Fq::from(3),
            "Wrong circuit evaluation"
        );

        // Layer poly checking output gate 01
        assert_eq!(
            LayerPoly::new(&vec![Fq::from(0), Fq::from(1)], &layer_op_mle, &layer_w_mle,).evaluate(
                &vec![
                    (1, Fq::from(0)),
                    (2, Fq::from(1)),
                    (3, Fq::from(0)),
                    (4, Fq::from(0)),
                    (5, Fq::from(1)),
                    (6, Fq::from(1))
                ]
            ),
            Fq::from(12),
            "Wrong circuit evaluation"
        );

        // Layer poly checking output gate 10
        assert_eq!(
            LayerPoly::new(&vec![Fq::from(1), Fq::from(0)], &layer_op_mle, &layer_w_mle,).evaluate(
                &vec![
                    (1, Fq::from(1)),
                    (2, Fq::from(0)),
                    (3, Fq::from(0)),
                    (4, Fq::from(1)),
                    (5, Fq::from(0)),
                    (6, Fq::from(1))
                ]
            ),
            Fq::from(11),
            "Wrong circuit evaluation"
        );
    }

    #[test]
    pub fn create_evaluate_layer_poly_at_a_random_point() {
        let circuit = create_circuit();

        let inputs = vec![
            Fq::from(1),
            Fq::from(2),
            Fq::from(3),
            Fq::from(4),
            Fq::from(5),
            Fq::from(6),
        ];

        let circuit_evaluation = circuit.evaluate(&inputs);

        // Create layer 0 poly
        let layer_op_mle = circuit.get_layer_ops_mles(0);

        // Layer w_mle for Layer 1 evaluations
        let layer_w_mle = get_wmle(&circuit_evaluation[1]);

        // w_mle for layer 0
        let layer_0_wmle = get_wmle(&circuit_evaluation[0]);

        // Layer poly checking output gate 5
        let layer_poly = LayerPoly::new(&vec![Fq::from(5)], &layer_op_mle, &layer_w_mle);

        assert_eq!(
            layer_poly.sum_over_the_boolean_hypercube(),
            layer_0_wmle.evaluate(&vec![(1, Fq::from(5))]),
            "Wrong layer poly evaluation"
        );

        // Layer 1
        let layer_op_mle = circuit.get_layer_ops_mles(1);

        // Layer w_mle for Layer 2 evaluations
        let layer_w_mle = get_wmle(&circuit_evaluation[2]);

        // w_mle for layer 1
        let layer_1_wmle = get_wmle(&circuit_evaluation[1]);

        // Layer poly checking output gate 10
        let layer_poly = LayerPoly::new(&vec![Fq::from(10)], &layer_op_mle, &layer_w_mle);

        assert_eq!(
            layer_poly.sum_over_the_boolean_hypercube(),
            layer_1_wmle.evaluate(&vec![(1, Fq::from(10))]),
            "Wrong layer poly evaluation"
        );

        // Layer 2
        let layer_op_mle = circuit.get_layer_ops_mles(2);

        // Layer w_mle for Layer 3 evaluations
        let layer_w_mle = get_wmle(&circuit_evaluation[3]);

        // w_mle for layer 2
        let layer_2_wmle = get_wmle(&circuit_evaluation[2]);

        // Layer poly checking output gate 10, 5
        let layer_poly = LayerPoly::new(
            &vec![Fq::from(10), Fq::from(5)],
            &layer_op_mle,
            &layer_w_mle,
        );

        assert_eq!(
            layer_poly.sum_over_the_boolean_hypercube(),
            layer_2_wmle.evaluate(&vec![(1, Fq::from(10)), (2, Fq::from(5))]),
            "Wrong layer poly evaluation"
        );
    }

    #[test]
    pub fn test_evaluate_layer_poly_at_a_random_point_with_alpha_and_beta() {
        let circuit = create_circuit_1();

        let inputs = vec![
            Fq::from(1),
            Fq::from(2),
            Fq::from(3),
            Fq::from(4),
            Fq::from(5),
            Fq::from(6),
        ];

        let circuit_evaluation = circuit.evaluate(&inputs);

        // Layer 2
        let layer_op_mle = circuit.get_layer_ops_mles(2);

        // Layer w_mle for Layer 3 evaluations
        let layer_w_mle = get_wmle(&circuit_evaluation[3]);

        // w_mle for layer 2
        let layer_2_wmle = get_wmle(&circuit_evaluation[2]);

        // Layer poly checking output gate 10, 5
        let layer_poly = LayerPoly::new(
            &vec![Fq::from(10), Fq::from(5)],
            &layer_op_mle,
            &layer_w_mle,
        );

        assert_eq!(
            layer_poly.sum_over_the_boolean_hypercube(),
            layer_2_wmle.evaluate(&vec![(1, Fq::from(10)), (2, Fq::from(5))]),
            "Wrong layer poly evaluation"
        );

        // using alpha = 3 and beta = 5
        let alpha = Fq::from(3);
        let beta = Fq::from(5);

        // using r_b = 10, 5 and r_c = 2,4
        let r_b = vec![(1, Fq::from(10)), (2, Fq::from(5))];
        let r_c = vec![(1, Fq::from(2)), (2, Fq::from(4))];

        let alpha_mul_wb_plus_beta_mul_wc =
            (alpha * layer_2_wmle.evaluate(&r_b)) + (beta * layer_2_wmle.evaluate(&r_c));

        let new_layer_op_mle: Vec<UniversalMLE<Fq>> = get_new_ops_mles(
            &layer_op_mle,
            &vec![alpha, beta],
            &vec![Fq::from(10), Fq::from(5)],
            &vec![Fq::from(2), Fq::from(4)],
        );

        let new_layer_mle = LayerPoly::new(&vec![], &new_layer_op_mle, &layer_w_mle);

        assert_eq!(
            alpha_mul_wb_plus_beta_mul_wc,
            new_layer_mle.poly.sum_over_the_boolean_hypercube(),
            "Incorrect results: Results do not match"
        );
    }

    #[test]
    pub fn test_evaluate_layer_mle_given_inputs() {
        let circuit = create_circuit();

        let layer_ops_mle = circuit.get_layer_ops_mles(2);

        let challenges = vec![
            Fq::from(0),
            Fq::from(1),
            Fq::from(0),
            Fq::from(1),
            Fq::from(0),
            Fq::from(0),
            Fq::from(1),
            Fq::from(1),
        ];

        let wb = Fq::from(3);

        let wc = Fq::from(4);

        let evaluation = evaluate_layer_mle_given_inputs(&layer_ops_mle, &challenges, &wb, &wc);

        assert_eq!(evaluation, Fq::from(12), "Evaluation does not match");
    }
}
