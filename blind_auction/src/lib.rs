use keccak_rust::*;

pub struct BlindAuction {
    pub commitments: Vec<[u8; 32]>,
}

pub trait Auction<T> {
    fn new() -> Self;
    fn commit(&mut self, commitment: &[u8; 32]) -> usize;
    fn open(&self, amount: usize, salt: &[u8], index: usize) -> Result<usize, String>;
}

impl Auction<usize> for BlindAuction {
    fn new() -> Self {
        BlindAuction {
            commitments: vec![],
        }
    }

    // Returns a commmitment and its index
    fn commit(&mut self, commitment: &[u8; 32]) -> usize {
        self.commitments.push(*commitment);
        self.commitments.len() - 1
    }

    fn open(&self, amount: usize, salt: &[u8], index: usize) -> Result<usize, String> {
        let mut hasher = Keccak::new(SecurityLevel::SHA256, StateBitsWidth::F1600);
        hasher.append(amount.to_be_bytes().as_slice());
        hasher.append(salt);
        let commitment_check = hasher.hash();
        let commitment = match self.commitments.get(index) {
            Some(val) => val,
            None => panic!("Commitment not found"),
        };
        assert_eq!(commitment_check, commitment, "Invalid commitment");
        Ok(amount)
    }
}

pub fn generate_commitment(amount: usize, salt: &[u8]) -> [u8; 32] {
    let mut hasher = Keccak::new(SecurityLevel::SHA256, StateBitsWidth::F1600);
    hasher.append(amount.to_be_bytes().as_slice());
    hasher.append(salt);
    hasher.hash().try_into().unwrap()
}

#[cfg(test)]
pub mod test {
    use crate::{generate_commitment, Auction, BlindAuction};

    #[test]
    fn test_blind_auction() {
        let mut auction = BlindAuction::new();

        let uche_commitment = generate_commitment(50, "Hello world".as_bytes());
        let uche_commitment_index = auction.commit(&uche_commitment);

        let abiodun_commitment = generate_commitment(100, "Hi there".as_bytes());
        let abiodun_commitment_index = auction.commit(&abiodun_commitment);

        let uche_reveal = auction
            .open(50, "Hello world".as_bytes(), uche_commitment_index)
            .unwrap();
        println!("Uche committed: {}", uche_reveal);

        let abiodun_reveal = auction
            .open(100, "Hi there".as_bytes(), abiodun_commitment_index)
            .unwrap();
        println!("Abiodun committed: {}", abiodun_reveal);
    }
}
