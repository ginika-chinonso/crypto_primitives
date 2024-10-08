# Fast Reed Solomon Interactive Oracle Proof of Proximity

The FRI protocol enables a prover prove to a verifier that a polynomial is of a low degree.
This library contains an implementation of the fri prover and verifier.

## STEPS
- Instantiate FRI protocol
```rust
let fri_protocol = Fri::new(KeccakHasher::new());
```

- Commit to a polynomial
```rust
let mut poly_commitment = fri_protocol.fri_commitment(poly, 32);
```

- FRI query phase
```rust
fri_protocol.fri_query(&mut poly_commitment, 32);
```

- FRI decommitment phase
```rust
let verify = fri_protocol.fri_decommitment(&poly_commitment, 32);
assert!(verify);
```

## Example
```rust
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

// Using a blowup factor of 2
let mut poly_commitment = fri_protocol.fri_commitment(trace_poly, trace.len() * 2);

fri_protocol.fri_query(
    &mut poly_commitment,
    trace.len() * 2,
);

let verify = fri_protocol.fri_decommitment(
    &poly_commitment,
    trace.len() * 2,
);

assert!(verify, "Invalid decommitment")
```