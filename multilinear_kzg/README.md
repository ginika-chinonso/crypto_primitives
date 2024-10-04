# Multilinear KZG Polynomial Commitment Scheme

This tool enables a prover commit to a multilinear polynomial and a verifier to open the polynomial at any point.
It also features a trusted set-up that can be contributed to

### STEPS
- Instantiate a trusted setup ceremony
```rust
let taus = vec![
    "contrib for a".as_bytes().to_vec(),
    "contrib for b".as_bytes().to_vec(),
];

let multilinear_pcs: MultilinearPCS<Bls12_381, DefaultFieldHasher<Sha256>> = MultilinearPCS::instantiate("initialize multilinear pcs".as_bytes(), taus);
```

- Contribute to the ceremony
```rust
coming soon!
```

- Commit to a polynomial
```rust
let poly = MLE::new(&vec![Fr::from(0), Fr::from(0), Fr::from(0), Fr::from(2)]);

let commitment = multilinear_pcs.commit(&poly);
```

- Open a polynomial
```rust
let proof = multilinear_pcs.open(&poly, &vec![(1, Fr::from(4)), (2, Fr::from(3))]);
```

- Verify proof
```rust
let verify = multilinear_pcs.verify(commitment, vec![Fr::from(4), Fr::from(3)], proof);
```

### Example
```rust
let taus = vec![
    "contrib for a".as_bytes().to_vec(),
    "contrib for b".as_bytes().to_vec(),
];

let multilinear_pcs: MultilinearPCS<Bls12_381, DefaultFieldHasher<Sha256>> = MultilinearPCS::instantiate("initialize multilinear pcs".as_bytes(), taus);

let poly = MLE::new(&vec![Fr::from(0), Fr::from(0), Fr::from(0), Fr::from(2)]);

let commitment = multilinear_pcs.commit(&poly);

let proof = multilinear_pcs.open(&poly, &vec![(1, Fr::from(4)), (2, Fr::from(3))]);

let verify = multilinear_pcs.verify(commitment, vec![Fr::from(4), Fr::from(3)], proof);

assert!(verify, "Invalid proof");
```