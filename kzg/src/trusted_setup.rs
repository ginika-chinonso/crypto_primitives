//////////////////////////////////////////
/// Trusted Setup Ceremony
//////////////////////////////////////////
use ark_ec::{pairing::Pairing, Group};
use ark_ff::{field_hashers::HashToField, Field, PrimeField};

#[derive(Debug)]
pub struct PublicParameters<C: Pairing> {
    // powers of tau for G1
    pub g1_powers_of_tau: Vec<C::G1>,
    // tau in G2
    pub g2_power_of_tau: Vec<C::G2>,
}

impl<C: Pairing> PublicParameters<C> {
    pub fn new(g1_powers_of_tau: Vec<C::G1>, g2_power_of_tau: Vec<C::G2>) -> Self {
        Self {
            g1_powers_of_tau,
            g2_power_of_tau,
        }
    }
}

pub struct TrustedSetUpCeremony<C: Pairing, H: HashToField<C::ScalarField>> {
    hasher: H,
    max_degree: u64,
    pub public_parameters: PublicParameters<C>,
    // Not considered for now
    // pub contributors_proof: Vec<F>
}

impl<C: Pairing, H: HashToField<C::ScalarField>> TrustedSetUpCeremony<C, H> {
    pub fn instantiate(hasher_domain: &[u8], degree: u64, randomness: &[u8]) -> Self {
        let hasher = HashToField::new(hasher_domain);

        let mut instance = Self {
            hasher,
            max_degree: degree,
            public_parameters: PublicParameters::new(vec![], vec![]),
        };

        let tau = instance
            .hasher
            .hash_to_field(randomness, 1)
            .get(0)
            .unwrap()
            .to_owned();

        let mut g1_powers_of_tau: Vec<C::G1> = vec![];

        let g = C::G1::generator();

        g1_powers_of_tau.push(g);

        for i in 0..degree {
            g1_powers_of_tau.push(g1_powers_of_tau[i as usize].mul_bigint(tau.into_bigint()));
        }

        let g2_power_of_tau = vec![C::G2::generator().mul_bigint(tau.pow([1]).into_bigint())];

        let public_parameters = PublicParameters::new(g1_powers_of_tau, g2_power_of_tau);

        instance.public_parameters = public_parameters;

        instance
    }

    pub fn contribute(&mut self, randomness: &[u8]) {
        let tau: C::ScalarField = self
            .hasher
            .hash_to_field(randomness, 1)
            .get(0)
            .unwrap()
            .to_owned();

        self.public_parameters.g1_powers_of_tau = self
            .public_parameters
            .g1_powers_of_tau
            .iter()
            .enumerate()
            .map(|(index, val)| val.mul_bigint(tau.pow([index as u64]).into_bigint()))
            .collect();

        self.public_parameters.g2_power_of_tau = self
            .public_parameters
            .g2_power_of_tau
            .iter()
            .map(|val| val.mul_bigint(tau.into_bigint()))
            .collect();
    }
}
