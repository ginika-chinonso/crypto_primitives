use ark_ec::{
    pairing::{Pairing, PairingOutput},
    Group,
};
use ark_ff::{field_hashers::HashToField, One, PrimeField, Zero};
use polynomials::multilinear_polynomial::{eval_form::MLE, traits::MultilinearPolynomialTrait};

pub struct MultilinearPCS<P: Pairing, H: HashToField<P::ScalarField>> {
    pub trusted_setup_parameter: TrustedSetup<P, H>,
}

#[derive(Debug, Clone)]
// Proof should be indexed by variable
pub struct MultilinearPCSProof<P: Pairing> {
    pub evaluation: P::G1,
    pub proof: Vec<P::G1>,
}

impl<P: Pairing> MultilinearPCSProof<P> {
    pub fn new(evaluation: P::G1, proof: Vec<P::G1>) -> Self {
        Self { evaluation, proof }
    }
}

impl<P: Pairing, H: HashToField<P::ScalarField>> MultilinearPCS<P, H> {
    pub fn instantiate(hasher_domain: &[u8], taus: Vec<Vec<u8>>) -> Self {
        let trusted_setup_parameter = TrustedSetup::new(hasher_domain, taus);

        MultilinearPCS {
            trusted_setup_parameter,
        }
    }

    pub fn contribute(&mut self, taus: Vec<Vec<u8>>) {
        self.trusted_setup_parameter.contribute(taus);

        self.trusted_setup_parameter.powers_of_tau =
            self.trusted_setup_parameter.powers_of_tau.clone();

        self.trusted_setup_parameter.tau_in_g2 = self.trusted_setup_parameter.tau_in_g2.clone();
    }

    pub fn commit(&self, poly: MLE<P::ScalarField>) -> P::G1 {
        self.trusted_setup_parameter
            .powers_of_tau
            .iter()
            .zip(poly.val)
            .map(|(tau, val)| tau.mul_bigint(val.into_bigint()))
            .fold(P::G1::zero(), |mut init, val| {
                init += val;
                init
            })
    }

    // Returns a vector of variables and its proof
    pub fn open(
        &self,
        poly: MLE<P::ScalarField>,
        eval_points: Vec<(usize, P::ScalarField)>,
    ) -> MultilinearPCSProof<P> {
        let evaluation = poly.evaluate(&eval_points);

        dbg!(&evaluation);

        let evaluation_in_g1 = P::G1::generator().mul_bigint(evaluation.into_bigint());

        let proof_polys = self.generate_multilinear_pcs_proof(
            poly.clone() - MLE::new(&vec![evaluation; poly.val.len()]),
            eval_points,
        );
        // let proof_polys = self.generate_multilinear_pcs_proof(poly.clone(), eval_points);

        let mut proof = vec![];

        // Commit to the proof polys
        for proof_poly in proof_polys.0 {
            dbg!(&proof_poly);
            proof.push(self.commit(proof_poly));
        }

        dbg!(&proof_polys.1);

        proof.push(P::G1::generator().mul_bigint(proof_polys.1.into_bigint()));

        dbg!(proof.len());

        MultilinearPCSProof::new(evaluation_in_g1, proof)
    }

    // Evaluation points should be indexed by variable
    pub fn verify(
        &self,
        commitment: P::G1,
        evaluation_points: Vec<P::ScalarField>,
        proof: MultilinearPCSProof<P>,
    ) -> bool {
        let lhs = commitment - proof.evaluation;

        dbg!(&evaluation_points.len());
        dbg!(&proof.proof.len());

        let rhs = evaluation_points
            .iter()
            .enumerate()
            .map(|(index, val)| {
                P::pairing(
                    proof.proof[index],
                    self.trusted_setup_parameter.tau_in_g2[index]
                        - P::G2::generator().mul_bigint(val.into_bigint()),
                )
            })
            .collect::<Vec<PairingOutput<P>>>();

        let final_sum_of_rhs = rhs.iter().skip(1).fold(
            rhs.get(0).expect("No proof available").clone(),
            |mut init, val| {
                init = init + val;
                init
            },
        );

        P::pairing(lhs, P::G2::generator()) == final_sum_of_rhs
    }

    pub fn generate_multilinear_pcs_proof(
        &self,
        mut poly: MLE<P::ScalarField>,
        mut evaluation_points: Vec<(usize, P::ScalarField)>,
    ) -> (Vec<MLE<P::ScalarField>>, P::ScalarField) {
        evaluation_points.sort();

        let mut remainders = (vec![], P::ScalarField::zero());

        // calculate the quotients
        let mut quotients = (vec![], P::ScalarField::zero());

        for i in 0..evaluation_points.len() {
            let quotient = poly.partial_eval(&vec![(1, P::ScalarField::one())])
                - poly.partial_eval(&vec![(1, P::ScalarField::zero())]);
            quotients
                .0
                .push(quotient.add_variable_at_index(&mut (1..=(i + 1)).collect::<Vec<usize>>()));
            poly = poly.partial_eval(&vec![(1, evaluation_points[i].1)]);
            remainders.0.push(poly.clone());
        }

        dbg!(&remainders);
        dbg!(&quotients);

        quotients
    }
}

////////////////////////////////////////////
/// Multilinear PCS Trusted Setup
/// ////////////////////////////////////////

#[derive(Debug, Clone)]
pub struct TrustedSetup<P: Pairing, H: HashToField<P::ScalarField>> {
    pub hasher: H,
    pub number_of_vars: usize,
    pub powers_of_tau: Vec<P::G1>,
    // Taus in G2 should be indexed by variable
    pub tau_in_g2: Vec<P::G2>,
}

impl<P, H> TrustedSetup<P, H>
where
    P: Pairing,
    H: HashToField<P::ScalarField>,
{
    pub fn new(hasher_domain: &[u8], taus: Vec<Vec<u8>>) -> Self {
        let hasher: H = HashToField::new(hasher_domain);

        let taus: Vec<P::ScalarField> = taus
            .iter()
            .map(|tau| {
                *hasher
                    .hash_to_field(tau, 1)
                    .get(0)
                    .expect("Error encountered when hashing to field")
            })
            .collect();

        let eq_for_tau_in_g1 = MLE::eq(&taus)
            .val
            .iter()
            .map(|val| P::G1::generator().mul_bigint(val.into_bigint()))
            .collect();

        let taus_in_g2 = taus
            .iter()
            .map(|val| P::G2::generator().mul_bigint(val.into_bigint()))
            .collect();

        Self {
            hasher,
            number_of_vars: taus.len(),
            powers_of_tau: eq_for_tau_in_g1,
            tau_in_g2: taus_in_g2,
        }
    }

    pub fn contribute(&mut self, randomness: Vec<Vec<u8>>) {
        assert_eq!(
            self.number_of_vars,
            randomness.len(),
            "Not enough taus for commitment"
        );

        let taus: Vec<P::ScalarField> = randomness
            .iter()
            .map(|val| {
                *self
                    .hasher
                    .hash_to_field(val, 1)
                    .get(0)
                    .expect("Error encountered when hashing to field")
            })
            .collect();

        let eq_for_taus = MLE::eq(&taus);

        self.powers_of_tau = self
            .powers_of_tau
            .iter()
            .zip(&eq_for_taus.val)
            .map(|(tau, val)| tau.mul_bigint(val.into_bigint()))
            .collect::<Vec<P::G1>>();

        self.tau_in_g2 = self
            .tau_in_g2
            .iter()
            .zip(&taus)
            .map(|(tau, val)| tau.mul_bigint(val.into_bigint()))
            .collect::<Vec<P::G2>>();
    }
}

//////////////////////////////
/// TEST
/// //////////////////////////

#[cfg(test)]
mod tests {
    use ark_bls12_381::{Bls12_381, Fr};
    use ark_ff::field_hashers::DefaultFieldHasher;
    use polynomials::multilinear_polynomial::eval_form::MLE;
    use sha2::Sha256;

    use crate::{MultilinearPCS, TrustedSetup};

    #[test]
    fn test_trusted_setup() {
        let taus: Vec<Vec<u8>> = vec![
            "tau for a".as_bytes().to_vec(),
            "tau for b".as_bytes().to_vec(),
        ];

        let mut trusted_setup: TrustedSetup<Bls12_381, DefaultFieldHasher<Sha256>> =
            TrustedSetup::new("hasher domain".as_bytes(), taus);

        trusted_setup.contribute(vec![
            "contribute to a".as_bytes().to_vec(),
            "contribute to b".as_bytes().to_vec(),
        ]);

        trusted_setup.contribute(vec![
            "second contribution to a".as_bytes().to_vec(),
            "second contribution to b".as_bytes().to_vec(),
        ]);

        dbg!(trusted_setup.powers_of_tau);
        dbg!(trusted_setup.tau_in_g2);
    }

    #[test]
    fn test_multilinear_pcs() {
        let taus = vec![
            "contrib for a".as_bytes().to_vec(),
            "contrib for b".as_bytes().to_vec(),
        ];

        let multilinear_pcs: MultilinearPCS<Bls12_381, DefaultFieldHasher<Sha256>> =
            MultilinearPCS::instantiate("initialize multilinear pcs".as_bytes(), taus);

        // TODO: Trouble shoot the contribute function
        // multilinear_pcs.contribute(vec!["second contrib for a".as_bytes().to_vec(), "second contrib for b".as_bytes().to_vec()]);

        // Polynomial in consideration: 2ab
        let val = vec![Fr::from(0), Fr::from(0), Fr::from(0), Fr::from(2)];

        let poly = MLE::new(&val);

        let commitment = multilinear_pcs.commit(poly.clone());

        dbg!(&commitment);

        let proof = multilinear_pcs.open(poly, vec![(1, Fr::from(4)), (2, Fr::from(3))]);

        dbg!(&proof);

        let verify = multilinear_pcs.verify(commitment, vec![Fr::from(4), Fr::from(3)], proof);

        dbg!(&verify);

        assert!(verify, "Invalid proof");
    }
}
