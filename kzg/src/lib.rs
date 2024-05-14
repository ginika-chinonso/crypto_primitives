use ark_ec::{pairing::Pairing, Group};
use ark_ff::{field_hashers::HashToField, fields::Field, PrimeField};

use polynomials::univariate_polynomial::UnivariatePolynomial;

// For a detailed explanation on KZG, visit: https://blog.subspace.network/kzg-polynomial-commitments-cd64af8ec868

// KZG

// if f(z) = b
// f(z) - b = 0
// Based on this, z is a root of this polynomial (f(z) - b)
// Therefore (x - z) is a factor of the polynomial f(z) - b
// Abstracting the value of z, we have
// f(x) - b
// where (x -z) is a factor
// Therefore, dividing f(x) - b by (x - z) should give us a polynomial without a remainder
// f(x) - b / (x - z) = q(x)
// where:
// f(x) = the polynomial we want to commit to
// x = can be any point we want to evaluate the polynomial at (in kzg it is used for the commitment and gotten from the trusted setup)
// z = the point we want to open the polynomial at
// b = the result of opening the polynomial at z
// q = quotient polynoial (should have no remainder)

pub struct KZG<C: Pairing, H: HashToField<C::ScalarField>> {
    trusted_setup: TrustedSetUpCeremony<C, H>,
}

// pub struct PolyCommitment<C: PrimeField > {
//     eval: C,
//     proof: C,
//     quotient_poly: UnivariatePolynomial<C>
// }

// impl <C: PrimeField > PolyCommitment<C> {

//     pub fn new(eval: C, proof: C, quotient_poly: UnivariatePolynomial<C>) -> Self {
//         Self { eval, proof, quotient_poly }
//     }
// }

// TODO:
// pub fn generate_zerifier_for_nth_roots_of_unity<F: PrimeField>(domain_size: usize) -> UnivariatePolynomial<F> {
//     // todo!()
// }
// pub fn generate_zerifier_for_domain<F: PrimeField>(domain_size: usize) -> UnivariatePolynomial<F> {
//     // todo!()
// }

impl<C: Pairing, H: HashToField<C::ScalarField>> KZG<C, H> {
    pub fn instantiate(
        hasher_domain: &[u8],
        domain_size: usize,
        initial_randomness: &[u8],
    ) -> Self {
        let trusted_setup_ceremony = TrustedSetUpCeremony::instantiate(
            hasher_domain,
            domain_size as u64,
            initial_randomness,
        );

        Self {
            trusted_setup: trusted_setup_ceremony,
        }
    }

    pub fn contribute(&mut self, randomness: &[u8]) {
        self.trusted_setup.contribute(randomness);
    }

    pub fn commit_to_poly(&self, poly: &UnivariatePolynomial<C::ScalarField>) -> C::G1 {
        assert!(
            self.trusted_setup.public_parameters.g1_powers_of_tau.len() >= poly.coefficients.len(),
            "Powers of tau not sufficient to evaluate poly"
        );

        let res = poly.coefficients.iter().enumerate().fold(
            C::G1::default(),
            |mut acc, (index, coeff)| {
                let res = self.trusted_setup.public_parameters.g1_powers_of_tau[index]
                    .mul_bigint(coeff.into_bigint());
                acc = acc + res;
                acc
            },
        );

        res
    }

    // Opens a polynomial (poly) at a point (eval_point) using the commitment to the polynomial
    // returns a tuple of the evaluation and proof
    pub fn open(
        &self,
        eval_point: C::ScalarField,
        poly: UnivariatePolynomial<C::ScalarField>,
        zeroifier: UnivariatePolynomial<C::ScalarField>,
    ) -> (C::ScalarField, C::G1) {
        assert!(
            self.trusted_setup.public_parameters.g1_powers_of_tau.len() >= poly.coefficients.len(),
            "Powers of tau not sufficient to evaluate poly"
        );

        let res = poly.evaluate(eval_point);

        let (quotient, _) = poly / zeroifier;

        let proof: C::G1 = quotient.coefficients.iter().enumerate().fold(
            C::G1::default(),
            |mut acc, (index, coeff)| {
                let res = self.trusted_setup.public_parameters.g1_powers_of_tau[index]
                    .mul_bigint(coeff.into_bigint());
                acc = acc + res;
                acc
            },
        );

        (res, proof)
    }

    pub fn verify(
        &self,
        eval_point: C::ScalarField,
        eval: C::ScalarField,
        commitment: C::G1,
        proof: C::G1,
    ) -> bool {
        let num = commitment - C::G1::generator().mul_bigint(eval.into_bigint());

        let denum = self.trusted_setup.public_parameters.g2_power_of_tau[0]
            - C::G2::generator().mul_bigint(eval_point.into_bigint());

        C::pairing(num, C::G2::generator()) == C::pairing(proof, denum)
    }
}

//////////////////////////////////////////
/// Trusted Setup Ceremony
//////////////////////////////////////////

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
    pub hasher: H,
    pub domain_size: u64,
    pub public_parameters: PublicParameters<C>, // Not considered for now
                                                // pub contributors_proof: Vec<F>
}

// TODO: Review
impl<C: Pairing, H: HashToField<C::ScalarField>> TrustedSetUpCeremony<C, H> {
    // TODO:
    // Correct domain size to degree
    pub fn instantiate(hasher_domain: &[u8], domain_size: u64, randomness: &[u8]) -> Self {
        let hasher = HashToField::new(hasher_domain);

        let mut instance = Self {
            hasher,
            domain_size,
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

        for i in 0..domain_size {
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

#[cfg(test)]
pub mod test {
    use crate::KZG;

    use super::TrustedSetUpCeremony;
    use ark_bls12_381::{Bls12_381, Fr};
    use ark_ec::CurveGroup;
    use ark_ff::field_hashers::DefaultFieldHasher;
    use polynomials::univariate_polynomial::UnivariatePolynomial;
    use sha2::Sha256;

    #[test]
    pub fn test_instantiate_trusted_setup() {
        let mut trusted_setup_instance: TrustedSetUpCeremony<
            Bls12_381,
            DefaultFieldHasher<Sha256>,
        > = TrustedSetUpCeremony::instantiate(
            "Hello world".as_bytes(),
            8,
            "Initialize Ceremony".as_bytes(),
        );

        trusted_setup_instance.contribute("My contribution".as_bytes());

        dbg!(trusted_setup_instance.public_parameters);
    }

    // TODO
    // Write a proper test
    #[test]
    pub fn test_commit_to_poly() {
        let poly = UnivariatePolynomial::new(vec![Fr::from(1), Fr::from(3), Fr::from(2)]);

        let kzg: KZG<Bls12_381, DefaultFieldHasher<Sha256>> = KZG::instantiate(
            "Hello world".as_bytes(),
            2,
            "Initialize Ceremony".as_bytes(),
        );

        let commitment = kzg.commit_to_poly(&poly).into_affine().to_string();

        dbg!(commitment);
    }

    #[test]
    pub fn test_kzg() {
        let poly = UnivariatePolynomial::new(vec![Fr::from(6), Fr::from(5), Fr::from(1)]);

        let mut kzg: KZG<Bls12_381, DefaultFieldHasher<Sha256>> = KZG::instantiate(
            "Hello world".as_bytes(),
            2,
            "Instantiate Ceremony".as_bytes(),
        );

        kzg.contribute("My contribution".as_bytes());

        kzg.contribute("Second contribution".as_bytes());

        let commitment = kzg.commit_to_poly(&poly);

        let zeroifier = UnivariatePolynomial::new(vec![-Fr::from(5000), Fr::from(1)]);

        let (evaluation, proof) = kzg.open(Fr::from(5000), poly, zeroifier);

        assert!(
            kzg.verify(Fr::from(5000), evaluation, commitment, proof),
            "Invalid proof"
        );
    }
}
