use ark_ec::{pairing::Pairing, Group};
use ark_ff::{field_hashers::HashToField, PrimeField};

use polynomials::univariate_polynomial::UnivariatePolynomial;

mod trusted_setup;
use trusted_setup::TrustedSetUpCeremony;

pub struct KZG<C: Pairing, H: HashToField<C::ScalarField>> {
    trusted_setup: TrustedSetUpCeremony<C, H>,
}

impl<C: Pairing, H: HashToField<C::ScalarField>> KZG<C, H> {
    fn instantiate(hasher_domain: &[u8], domain_size: usize, initial_randomness: &[u8]) -> Self {
        let trusted_setup_ceremony = TrustedSetUpCeremony::instantiate(
            hasher_domain,
            domain_size as u64,
            initial_randomness,
        );

        Self {
            trusted_setup: trusted_setup_ceremony,
        }
    }

    fn contribute(&mut self, randomness: &[u8]) {
        self.trusted_setup.contribute(randomness);
    }

    fn commit_to_poly(&self, poly: &UnivariatePolynomial<C::ScalarField>) -> C::G1 {
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
    fn open(
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

    fn verify(
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
