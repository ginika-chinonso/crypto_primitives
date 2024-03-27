use keccak_rust::*;

#[derive(Debug, Clone)]
pub struct MerkleTree {
    pub tree: Vec<Vec<[u8; 32]>>,
}

pub trait Hasher {
    fn new() -> Self;
    fn hash(&mut self, input: &[u8]) -> [u8; 32];
    fn hash_two_leaves(&mut self, leaf1: &[u8; 32], leaf2: &[u8; 32]) -> [u8; 32];
}

// struct PoseidonHasher(Poseidon);

// impl Hasher for PoseidonHasher {
//     fn new() -> Self {
//         PoseidonHasher(Poseidon::new())
//     }

//     fn hash<F: PrimeField>(&self, input: &[u8]) -> [u8; 32] {
//         let input_to_field_element = F::from_random_bytes(input).unwrap();
//         let res: [u8; 32] = self.0.hash(vec![input_to_field_element]).unwrap();
//         res
//     }

// }

struct KeccakHasher(Keccak);

impl Hasher for KeccakHasher {
    fn new() -> Self {
        KeccakHasher(Keccak::new(SecurityLevel::SHA256, StateBitsWidth::F1600))
        // Keccak::new(keccak_rust::SecurityLevel::SHA256, StateBitsWidth::F1600))
    }

    fn hash(&mut self, input: &[u8]) -> [u8; 32] {
        self.0.append(input);
        self.0.hash().try_into().unwrap()
    }

    fn hash_two_leaves(&mut self, leaf1: &[u8; 32], leaf2: &[u8; 32]) -> [u8; 32] {
        self.hash(leaf1);
        self.hash(leaf2)
    }
}

impl MerkleTree {
    pub fn new<H: Hasher>(leaf_node: &Vec<Vec<u8>>) -> Self {
        let mut res = vec![];

        let hashed_leaves: Vec<[u8; 32]> = leaf_node
            .into_iter()
            .map(|leaf| H::new().hash(leaf))
            .collect();

        let leaves = if hashed_leaves.len() % 2 == 0 {
            hashed_leaves.to_vec()
        } else {
            // we need to pad the leaves to make it even
            // we use the last element in the leaves to do the padding
            let mut new_leaves = hashed_leaves.to_vec();
            new_leaves.push(hashed_leaves[hashed_leaves.len() - 1]);
            new_leaves
        };

        res.push(leaves);

        loop {
            if res[res.len() - 1].len() == 2 {
                let mut hasher = H::new();
                res.push(vec![
                    hasher.hash_two_leaves(&res[res.len() - 1][0], &res[res.len() - 1][1])
                ]);
                break;
            }

            let mut next_layer = vec![];

            for i in (0..res[res.len() - 1].len()).step_by(2) {
                let mut hasher = H::new();

                let parent_hash =
                    hasher.hash_two_leaves(&res[res.len() - 1][i], &res[res.len() - 1][i + 1]);

                next_layer.push(parent_hash);
            }

            if next_layer.len() % 2 == 0 {
                res.push(next_layer);
            } else {
                next_layer.push(next_layer[next_layer.len() - 1]);
                res.push(next_layer);
            }
        }

        res.reverse();

        Self { tree: res }
    }

    pub fn get_root(&self) -> Option<[u8; 32]> {
        match self.tree[0].get(0) {
            Some(root) => Some(*root),
            None => None,
        }
    }

    // Should this take the index or the leaf?
    pub fn get_proof(&self, index: usize) -> Vec<[u8; 32]> {
        let mut res: Vec<[u8; 32]> = vec![];
        let tree_depth= self.tree_depth();
        let mut current_index = index;
        res.push(self.tree[tree_depth - 1][get_sibling_index(current_index)]);

        // this is incorrect. rectify
        for i in (1..tree_depth - 1).rev() {

            let parent_index = get_parent_index(current_index);
            let parent_sibling = get_sibling_index(parent_index);
            res.push(self.tree[i][parent_sibling]);
            // res.push(self.tree[i - 1][get_parent_index(current_index)]);
            current_index = parent_index;
        }
        res.reverse();
        res
    }

    fn tree_depth(&self) -> usize {
        self.tree.len()
    }
}

// Verify
pub fn verify_proof<H: Hasher, T>(
    root_hash: [u8; 32],
    index: usize,
    leaf: &[u8],
    proof: Vec<[u8; 32]>,
) -> bool {
    let mut leaf_hash = H::new().hash(leaf);
    let mut current_index = index;

    for i in (0..proof.len()).rev() {
        let mut hasher = H::new();
        if current_index % 2 == 0 {
            dbg!(&current_index);
            leaf_hash = hasher.hash_two_leaves(&leaf_hash, &proof[i]);
            current_index = get_parent_index(current_index);
            dbg!(&current_index);
        } else {
            dbg!(&current_index);
            leaf_hash = hasher.hash_two_leaves(&proof[i], &leaf_hash);
            current_index = get_parent_index(current_index);
            dbg!(&current_index);
        };
    }
    root_hash == leaf_hash
}

fn get_parent_index(index: usize) -> usize {
    index / 2
}

fn get_sibling_index(index: usize) -> usize {

    if index % 2 == 0 {
        index + 1
    } else {
        index - 1
    }
}






#[cfg(test)]
pub mod test {

    use crate::{verify_proof, Hasher, KeccakHasher, MerkleTree};

    #[test]
    pub fn test_create_merkle_tree() {
        let array: Vec<Vec<u8>> = vec![1, 2, 3]
            .into_iter()
            .map(|val: i32| val.to_be_bytes().to_vec())
            .collect();

        let merkle_tree = MerkleTree::new::<KeccakHasher>(&array);

        assert!(
            merkle_tree.tree[2][0] == KeccakHasher::new().hash(&1_i32.to_be_bytes().to_vec()),
            "Invalid leaf"
        );
        assert!(
            merkle_tree.tree[2][1] == KeccakHasher::new().hash(&2_i32.to_be_bytes().to_vec()),
            "Invalid leaf"
        );
        assert!(
            merkle_tree.tree[2][2] == KeccakHasher::new().hash(&3_i32.to_be_bytes().to_vec()),
            "Invalid leaf"
        );
        assert!(
            merkle_tree.tree[2][3] == KeccakHasher::new().hash(&3_i32.to_be_bytes().to_vec()),
            "Invalid leaf"
        );

        let one_hash = KeccakHasher::new().hash(1_i32.to_be_bytes().as_slice());
        let two_hash = KeccakHasher::new().hash(2_i32.to_be_bytes().as_slice());
        let three_hash = KeccakHasher::new().hash(3_i32.to_be_bytes().as_slice());

        assert!(
            merkle_tree.tree[1][0] == KeccakHasher::new().hash_two_leaves(&one_hash, &two_hash),
            "Invalid node"
        );
        assert!(
            merkle_tree.tree[1][1] == KeccakHasher::new().hash_two_leaves(&three_hash, &three_hash),
            "Invalid node"
        );

        let root_hash = KeccakHasher::new().hash_two_leaves(
            &KeccakHasher::new().hash_two_leaves(&one_hash, &two_hash),
            &KeccakHasher::new().hash_two_leaves(&three_hash, &three_hash),
        );

        assert!(merkle_tree.tree[0][0] == root_hash, "Invalid root");
    }

    #[test]
    pub fn test_generate_proof() {
        let array: Vec<Vec<u8>> = vec![1, 2, 3]
            .into_iter()
            .map(|val: i32| val.to_be_bytes().to_vec())
            .collect();

        let merkle_tree = MerkleTree::new::<KeccakHasher>(&array);

        let proof = merkle_tree.get_proof(2);

        assert!(proof.len() == 2, "Invalid proof length");

        assert!(proof[proof.len() - 1] == KeccakHasher::new().hash(3_i32.to_be_bytes().as_slice()));

        let one_hash = KeccakHasher::new().hash(1_i32.to_be_bytes().as_slice());
        let two_hash = KeccakHasher::new().hash(2_i32.to_be_bytes().as_slice());

        assert!(proof[proof.len() - 2] == KeccakHasher::new().hash_two_leaves(&one_hash, &two_hash));

    }

    #[test]
    pub fn test_proof_and_verify() {
        let array = (1..300).into_iter().map(|val| val)
            .into_iter()
            .map(|val: i32| val.to_be_bytes().to_vec())
            .collect();

        let merkle_tree = MerkleTree::new::<KeccakHasher>(&array);

        let proof = merkle_tree.get_proof(2);

        assert!(verify_proof::<KeccakHasher, i32>(
            merkle_tree.get_root().unwrap(),
            2,
            3_i32.to_be_bytes().as_slice(),
            proof
        ), "Invalid proof");
    }
}
