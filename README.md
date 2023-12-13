# Crypto Primitives
This repository is a random experiment at naively implementing some crypto primitives while learning rust

## Overview
Things to implement include;

### - [ ] Field
    - [v] New
    - [v] Modulus
### - [ ] Field Element
    - [v] New
    - [v] Value
    - [v] Field
    - [v] Add
    - [ ] Subtract
    - [ ] Multiply
    - [ ] Divide
    - [ ] Inverse
### - [ ] Elliptic curve
    - [v] New
    - [ ] Is on curve
    - [ ] Add
    - [ ] Multiply
### - [ ] Elliptic curve point
    - [v] New
    - [v] Is on curve
    - [v] Add
    - [v] Double point
    - [v] Scalar mul
### - [ ] Univariate Polynomial
    - [v] New
    - [v] Add
    - [v] Multiply
    - [v] Evaluate
    - [v] Interpolate
    - [v] Is zero
    - [v] Truncate
    - [v] From // converts from multilinear to univariate
    - [ ] Sparse to dense
### - [ ] Multilinear Polynomial
    - [v] New
    - [v] Truncate
    - [v] Simplify
    - [v] Add
    - [v] Partial evaluate
    - [v] Scalar mul
    - [v] Multiply
    - [v] Evaluate
    - [v] Interpolate
    - [v] Relabel
    - [ ] From // convert from univariate to multilinear
    - [ ] Is zero
    - [ ] Sparse to dense
### - [ ] Protocols
    - [v] Reed-Solomon Fingerprinting
    - [ ] Sumcheck protocol
    - [ ] GKR protocol
    - [ ] Fiat Shamir
    - [ ] Inner Product Argument (IPA)
    - [ ] KZG
    - [ ] Bullet Proofs
    - [ ] Groth 16

### - [ ] Benchmarking
    - [v] Polynomial ( Univariate )
    - [v] Multivariate Polynomial
