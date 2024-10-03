use ark_ec::{pairing::Pairing, Group};
use ark_ff::{field_hashers::HashToField, PrimeField};
use polynomials::multilinear_polynomial::eval_form::MLE;

////////////////////////////////////////////
/// Multilinear PCS Trusted Setup
/// ////////////////////////////////////////

#[derive(Debug, Clone)]
pub struct TrustedSetup<P: Pairing, H: HashToField<P::ScalarField>> {
    hasher: H,
    number_of_vars: usize,
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

    pub fn contribute(&mut self, randomness: &Vec<Vec<u8>>) {
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
