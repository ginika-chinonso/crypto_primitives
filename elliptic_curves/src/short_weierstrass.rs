use std::ops::{Add, BitXor};

use crate::fields::arithmetics::FieldElement;

// An elliptic curve is denoted by the equation
// y^2 = x^3 + ax + b -> Short weierstras
// y^2 + a1xy + a3y = x^3 + a2x^2 + a4x + a6 -> General weierstras
// If 4a^3 + 27b^2 != 0 the curve is non singular
#[derive(Debug, Clone, Copy)]
pub struct EllipticCurve {
    a: FieldElement,
    b: FieldElement,
}

impl EllipticCurve {
    pub fn new(a: FieldElement, b: FieldElement) -> EllipticCurve {
        EllipticCurve { a, b }
    }

    pub fn ec_point(&self, x: FieldElement, y: FieldElement) -> ECPoint {
        ECPoint { x, y, curve: *self }
    }
}

#[derive(Debug, Clone, Copy)]
pub struct ECPoint {
    x: FieldElement,
    y: FieldElement,
    curve: EllipticCurve,
}

impl ECPoint {
    pub fn new(x: FieldElement, y: FieldElement, curve: EllipticCurve) -> Self {
        Self { x, y, curve }
    }

    // pub fn is_on_curve(&self, curve: EllipticCurve) -> bool {
    //     if self.a. != curve.a || self.b != curve.b {
    //         false
    //     } else {
    //         true
    //     }
    // }

    // For point doubling
    // s = (3x1^2 + a) / 2y1
    // x3 = s^2 - 2x1
    // y3 = s(x1 - x3) - y1
    // pub fn double_point(self) -> Self {
    //     let num: F = (F::from(3) * self.x * self.x) + self.a;
    //     let denom: F = F::from(2) * self.y;
    //     let slope = num / denom;
    //     let x3 = (slope * slope) - (F::from(2) * self.x);
    //     let y3 = slope * (self.x - x3) - self.y;
    //     Self {
    //         x: x3,
    //         y: y3,
    //         curve: self.curve
    //     }
    // }

    // pub fn scalar_mul(self, scalar: usize) -> Self {
    //     for _i in 0..scalar {
    //         self.double_point();
    //     }
    //     self
    // }
}

// impl Add for ECPoint {
//     type Output = ECPoint;

// For point addition
// s = (y2 - y1) / (x2 - x1)
// x3 = s^2 - x1 - x2
// y3  = s(x1 - x3) - y1
// fn add(self, rhs: Self) -> Self::Output {

//     assert!(self.a == p2.a && self.b == p2.b);
//     let slope = (p2.y - self.y) / (p2.x - self.x);
//     let x3 = (slope * slope) - p2.x - self.x;
//     let y3 = slope * (p2.x - x3) - p2.y;
//     Self {
//         x: x3,
//         y: y3,
//         curve: self.curve
//     }
// }
// }

// impl BitXor for ECPoint {
//     type Output = ECPoint;

//     fn bitxor(self, rhs: usize) -> Self::Output {
//         todo!()
//     }
// }

#[cfg(test)]
mod tests {
    use super::EllipticCurve;
    use ark_ff::{Fp64, MontBackend, MontConfig};

    #[derive(MontConfig)]
    #[modulus = "17"]
    #[generator = "3"]
    pub struct FqConfig;
    pub type Fq = Fp64<MontBackend<FqConfig, 1>>;

    // #[test]
    // fn test_create_ec_point() {
    //     let _new_curve = EllipticCurve::new(Fq::from(7), Fq::from(20));
    // }
}
