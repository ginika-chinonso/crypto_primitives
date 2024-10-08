use itertools::{EitherOrBoth, Itertools};
use merkle_tree::{verify_proof, Hasher, KeccakHasher, MerkleTree};
use std::marker::PhantomData;

use ark_ff::{BigInteger, FftField, PrimeField};
use fiat_shamir_transcript::Transcript;
use polynomials::univariate_polynomial::UnivariatePolynomial;

#[derive(Debug)]
pub struct FriProof<F: PrimeField> {
    // The last round poly can as well be a constant depending on the implementation
    round_polys: Vec<UnivariatePolynomial<F>>,
    rounds_poly_commitment: Vec<[u8; 32]>,
    rounds_merkle_tree: Vec<MerkleTree>,
    query_indices: Vec<usize>,
    rounds_decommitment: Vec<FriDecommitment<F>>,
}

impl<F: PrimeField> FriProof<F> {
    pub fn new() -> Self {
        FriProof {
            round_polys: vec![],
            rounds_poly_commitment: vec![],
            rounds_decommitment: vec![],
            rounds_merkle_tree: vec![],
            query_indices: vec![],
        }
    }
}

#[derive(Debug)]
pub struct FriDecommitment<F: PrimeField> {
    value: F,
    sib: F,
    value_merkle_authentication_path: Vec<[u8; 32]>,
    sib_merkle_authentication_path: Vec<[u8; 32]>,
}

impl<F: PrimeField> FriDecommitment<F> {
    pub fn new(
        value: F,
        sib: F,
        value_merkle_authentication_path: Vec<[u8; 32]>,
        sib_merkle_authentication_path: Vec<[u8; 32]>,
    ) -> Self {
        Self {
            value,
            sib,
            value_merkle_authentication_path,
            sib_merkle_authentication_path,
        }
    }
}

pub struct Fri<F, H> {
    _marker: PhantomData<F>,
    hasher: H,
}

impl<F: FftField + PrimeField + std::convert::From<i32>, H: Hasher> Fri<F, H> {
    // Instantiate a FRI protocol
    pub fn new(hasher: H) -> Fri<F, H> {
        Fri {
            _marker: PhantomData,
            hasher,
        }
    }

    fn fold_poly(&self, beta: F, poly: &UnivariatePolynomial<F>) -> UnivariatePolynomial<F> {
        let (even_poly, odd_poly) =
            poly.coefficients
                .iter()
                .enumerate()
                .fold((vec![], vec![]), |mut acc, (index, val)| {
                    if index % 2 == 0 {
                        acc.0.push(val);
                        (acc.0, acc.1)
                    } else {
                        acc.1.push(val);
                        (acc.0, acc.1)
                    }
                });

        let reduced_odd_poly: Vec<F> = odd_poly.into_iter().map(|coef| *coef * beta).collect();

        let new_coefficients =
            even_poly
                .into_iter()
                .zip_longest(reduced_odd_poly)
                .fold(vec![], |mut acc, val| {
                    match val {
                        EitherOrBoth::Both(a, b) => acc.push(*a + b),
                        EitherOrBoth::Left(a) => acc.push(*a + F::zero()),
                        EitherOrBoth::Right(b) => acc.push(b + F::zero()),
                    };
                    acc
                });

        UnivariatePolynomial::new(new_coefficients)
    }

    fn poly_commitment(
        &self,
        poly: &UnivariatePolynomial<F>,
        commitment_domain: u64,
    ) -> MerkleTree {
        let eval_domain: F = generate_eval_domain(commitment_domain).unwrap();

        let mut evaluations = vec![F::zero(); commitment_domain as usize];

        for i in 0..commitment_domain {
            evaluations[i as usize] = poly.evaluate(eval_domain.pow([i]));
        }

        let leaf_node = evaluations
            .into_iter()
            .map(|val| val.into_bigint().to_bytes_be())
            .collect();

        let merkle_tree = MerkleTree::new::<H>(&leaf_node);

        merkle_tree
    }

    pub fn fri_commitment(&self, poly: UnivariatePolynomial<F>, trace_length: u64) -> FriProof<F> {
        let mut res = FriProof::new();

        let mut transcript = Transcript::new();

        let mut trace_length = trace_length;

        let mut round_poly = poly.clone();

        // find an better way of getting the number of rounds
        for _ in 0..((poly.coefficients.len() as f64).log2() as usize) {
            transcript.add_univariate_poly(&round_poly);

            // commit to the polynomial in evaluation form
            let round_commitment = self.poly_commitment(&round_poly, trace_length);
            res.round_polys.push(round_poly.clone());
            res.rounds_poly_commitment
                .push(round_commitment.get_root().unwrap());
            res.rounds_merkle_tree.push(round_commitment);

            let beta = transcript.sample_field_element();
            round_poly = self.fold_poly(beta, &round_poly);
            trace_length /= 2;
        }

        res
    }

    pub fn fri_query(&self, fri_proof: &mut FriProof<F>, domain_size: usize) {
        let mut domain_size = domain_size as u64;

        let transcript = Transcript::new();

        for i in 0..fri_proof.round_polys.len() {
            let index_to_query_at = transcript.sample_element_in_range_n(domain_size);

            let index_sib = (index_to_query_at + (domain_size >> 1)) % domain_size;

            let round_merkle_tree = &fri_proof.rounds_merkle_tree[i];

            let val_proof = round_merkle_tree.get_proof(index_to_query_at as usize);
            let sib_proof = round_merkle_tree.get_proof(index_sib as usize);

            let val = fri_proof.round_polys[i].evaluate(
                generate_eval_domain::<F>(domain_size)
                    .unwrap()
                    .pow([index_to_query_at]),
            );

            let sib = fri_proof.round_polys[i].evaluate(
                generate_eval_domain::<F>(domain_size)
                    .unwrap()
                    .pow([index_sib]),
            );

            let decommitment = FriDecommitment::new(val, sib, val_proof, sib_proof);

            fri_proof.rounds_decommitment.push(decommitment);
            fri_proof.query_indices.push(index_to_query_at as usize);

            domain_size >>= 1;
        }
    }

    pub fn fri_decommitment(&self, proof: &FriProof<F>, domain_size: usize) -> bool {
        let mut domain_size = domain_size;

        let mut transcript = Transcript::new();

        for i in 0..proof.rounds_decommitment.len() - 1 {
            transcript.add_univariate_poly(&proof.round_polys[i]);

            let beta: F = transcript.sample_field_element();

            let round_decommitment = &proof.rounds_decommitment[i];

            let cp_i_plus_one_x = &proof.round_polys[i + 1];

            let round_merkle_tree = &proof.rounds_merkle_tree[i];

            let index = proof.query_indices[i];

            let sib_index = proof.query_indices[i] + ((domain_size >> 1) % domain_size) as usize;

            let ex = generate_eval_domain::<F>(domain_size as u64)
                .unwrap()
                .pow([index as u64]);

            let minus_ex = generate_eval_domain::<F>(domain_size as u64)
                .unwrap()
                .pow([sib_index as u64]);

            assert!(
                verify_proof::<KeccakHasher>(
                    round_merkle_tree.get_root().unwrap(),
                    index,
                    &round_decommitment.value.into_bigint().to_bytes_be(),
                    round_decommitment.value_merkle_authentication_path.clone()
                ),
                "Merkle proof did not match"
            );
            assert!(
                verify_proof::<KeccakHasher>(
                    round_merkle_tree.get_root().unwrap(),
                    sib_index,
                    &round_decommitment.sib.into_bigint().to_bytes_be(),
                    round_decommitment.sib_merkle_authentication_path.clone()
                ),
                "Merkle proof did not match"
            );

            let g = ((beta + ex) / (ex.mul(F::from(2)))) * round_decommitment.value;

            let h = ((beta + minus_ex) / (minus_ex.mul(F::from(2)))) * round_decommitment.sib;

            assert_eq!(
                cp_i_plus_one_x.evaluate(ex.pow(F::from(2).into_bigint())),
                g + h,
                "Verifier check failed"
            );

            domain_size >>= 1;
        }

        true
    }
}

fn generate_eval_domain<F: FftField + PrimeField>(n: u64) -> Option<F> {
    F::get_root_of_unity(n)
}

#[cfg(test)]
pub mod test {
    use ark_ff::{BigInteger, Field, PrimeField};
    use merkle_tree::{Hasher, KeccakHasher, MerkleTree};
    use polynomials::univariate_polynomial::UnivariatePolynomial;

    use super::{generate_eval_domain, Fri};

    use ark_bn254::Fr;
    type Fq = Fr;

    // use ark_ff::{Fp64, MontBackend, MontConfig};
    // #[derive(MontConfig)]
    // #[modulus = "17"]
    // #[generator = "3"]
    // pub struct FqConfig;
    // pub type Fq = Fp64<MontBackend<FqConfig, 1>>;

    #[test]
    pub fn test_reduce_poly_degree() {
        let poly = UnivariatePolynomial::new(vec![
            Fq::from(6),
            Fq::from(2),
            Fq::from(3),
            Fq::from(9),
            Fq::from(1),
            Fq::from(1),
            Fq::from(1),
        ]);

        let fri: Fri<Fq, KeccakHasher> = Fri::new(KeccakHasher::new());

        let beta = Fq::from(5);

        let reduced_poly = fri.fold_poly(beta, &poly);

        // even_powers = [6, 3, 1, 1]
        // odd_powers = [2, 9, 1]
        // scaled_odd = [10, 45, 5]
        // new_poly = [16, 48, 6, 1]

        assert!(
            reduced_poly
                == UnivariatePolynomial::new(vec![
                    Fq::from(16),
                    Fq::from(48),
                    Fq::from(6),
                    Fq::from(1)
                ])
        );
    }

    #[test]
    pub fn test_generate_eval_domain() {
        let trace_length = 16_u64;
        let eval_domain = generate_eval_domain::<Fq>(trace_length).unwrap();
        let res = eval_domain.pow(&[trace_length]);
        assert!(res == Fq::from(1));
    }

    #[test]
    pub fn test_fri_commitment() {
        let poly = UnivariatePolynomial::new(vec![
            Fq::from(6),
            Fq::from(2),
            Fq::from(3),
            Fq::from(9),
            Fq::from(1),
            Fq::from(1),
            Fq::from(1),
        ]);

        let fri: Fri<Fq, KeccakHasher> = Fri::new(KeccakHasher::new());

        let poly_commitment: MerkleTree = fri.poly_commitment(&poly, 8);

        let root_of_unity_generator: Fq = generate_eval_domain(8).unwrap();

        let query_index = 3;

        let poly_eval_at_query: Fq = poly.evaluate(root_of_unity_generator.pow([query_index]));

        assert!(
            KeccakHasher::new().hash(poly_eval_at_query.into_bigint().to_bytes_be().as_slice())
                == poly_commitment.get_leaf_at_index(query_index),
            "Wrong evaluation"
        );
    }

    #[test]
    pub fn test_fri_protocol() {
        let trace = vec![
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
        ];

        let generator: Fq = generate_eval_domain(trace.len() as u64).unwrap();

        let domain: Vec<Fq> = (0..trace.len())
            .into_iter()
            .map(|index| generator.pow([index as u64]))
            .collect();

        let trace_poly = UnivariatePolynomial::interpolate(&domain, &trace);

        let fri_protocol = Fri::new(KeccakHasher::new());

        // // Using a blowup factor of 2
        let mut poly_commitment = fri_protocol.fri_commitment(trace_poly, (trace.len() * 2) as u64);

        fri_protocol.fri_query(&mut poly_commitment, trace.len() * 2);

        let decommitment = fri_protocol.fri_decommitment(&poly_commitment, trace.len() * 2);

        assert!(decommitment, "Invalid decommitment")
    }

    #[test]
    fn test_decommitment_algo_1() {
        let poly = UnivariatePolynomial::new(vec![
            Fq::from(1),
            Fq::from(2),
            Fq::from(3),
            Fq::from(2),
            Fq::from(3),
            Fq::from(4),
        ]);

        let fri: Fri<Fq, KeccakHasher> = Fri::new(KeccakHasher::new());

        let beta = Fq::from(3);

        let reduced_poly = fri.fold_poly(beta, &poly);

        assert!(
            reduced_poly == UnivariatePolynomial::new(vec![Fq::from(7), Fq::from(9), Fq::from(15)]),
            "Incorrect reduced poly"
        );

        let ex = Fq::from(9);
        let minus_ex = Fq::from(-9);

        let cp_i_x = poly;
        let cp_i_plus_one_x = reduced_poly;

        let g = ((beta + ex) / (ex * (Fq::from(2)))) * cp_i_x.evaluate(ex);

        let h = ((beta + minus_ex) / (minus_ex * (Fq::from(2)))) * cp_i_x.evaluate(minus_ex);

        assert_eq!(
            cp_i_plus_one_x.evaluate(ex.pow(Fq::from(2).into_bigint())),
            g + h,
            "Verifier check failed"
        );
    }

    #[test]
    fn test_decommitment_algo_2() {
        let poly = UnivariatePolynomial::new(vec![
            Fq::from(6),
            Fq::from(2),
            Fq::from(3),
            Fq::from(9),
            Fq::from(1),
            Fq::from(1),
            Fq::from(1),
        ]);

        let fri: Fri<Fq, KeccakHasher> = Fri::new(KeccakHasher::new());

        let beta = Fq::from(3);

        let reduced_poly = fri.fold_poly(beta, &poly);

        // even_powers = [6, 3, 1, 1]
        // odd_powers = [2, 9, 1]
        // scaled_odd = [6, 27, 3]
        // new_poly = [12, 30, 4, 1]

        assert!(
            reduced_poly
                == UnivariatePolynomial::new(vec![
                    Fq::from(12),
                    Fq::from(30),
                    Fq::from(4),
                    Fq::from(1)
                ]),
            "Incorrect reduced poly"
        );

        let ex = Fq::from(500);
        let minus_ex = Fq::from(-500);

        let g = ((beta + ex) / (ex * (Fq::from(2)))) * poly.evaluate(ex);

        let h = ((beta + minus_ex) / (minus_ex * (Fq::from(2)))) * poly.evaluate(minus_ex);

        assert_eq!(
            reduced_poly.evaluate(ex.pow(Fq::from(2).into_bigint())),
            g + h,
            "Verifier check failed"
        );
    }
}
