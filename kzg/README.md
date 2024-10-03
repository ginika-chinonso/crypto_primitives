# KZG Polynomial Commitment Scheme

This tool enables a prover commit to a polynomial and a verifier to open the polynomial at any point.
It also features a trusted set-up that can be contributed to

### STEPS
- Instantiate a trusted setup ceremony
```rust
let mut kzg: KZG<Bls12_381, DefaultFieldHasher<Sha256>> = KZG::instantiate(
    "Hello world".as_bytes(),
    2,
    "Instantiate Ceremony".as_bytes(),
)
```

- Contribute to the ceremony
```rust
kzg.contribute("My contribution".as_bytes());
```

- Commit to a polynomial
```rust
let commitment = kzg.commit_to_poly(&poly);
```

- Open a polynomial
```rust
let (evaluation, proof) = kzg.open(Fr::from(5000), poly, zeroifier);
```

- Verify proof
```rust
kzg.verify(Fr::from(5000), evaluation, commitment, proof),
```

### Example
```rust
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
```