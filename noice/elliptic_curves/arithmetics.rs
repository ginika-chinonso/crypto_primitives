
use crate::fields::arithmetics::{*};


// An elliptic curve is denoted by the equation
// y^2 = x^3 + ax + b
#[derive(Debug, Clone)]
pub struct EC<F, FE>{
    a: FE,
    b: FE,
    pub field: F,
}

#[derive(Debug)]
pub struct FieldPoint<FE> {
    x: FE,
    y: FE,
    a: FE,
    b: FE
}

pub trait ECTrait<F, FE>{
    fn new(a: FE, b: FE, field: F) -> EC<F, FE>;
    fn is_on_curve(self: Self, point: FieldPoint<FE>) -> bool;
    // fn add_points(p1: FieldPoint<FE>, p2: FieldPoint<FE>) -> FieldPoint<FE>;
    // fn scalar_mul(scalar: i32, point: FieldPoint<FE>) -> FieldPoint<FE>;
}



impl<F: FieldTrait + Clone + Copy, FE: FieldElementTrait<F>> ECTrait<F, FE> for EC<F, FE> {
    fn new(a: FE, b: FE, field: F) -> EC<F, FE> {
        EC { a, b, field }
    }

    fn is_on_curve(self: Self, point: FieldPoint<FE>) -> bool {
        if point.a.value() != self.a.value() || point.b.value() != self.b.value() {
            false
        } else {true}
    }

    // fn add_points(p1: FieldPoint<FE>, p2: FieldPoint<FE>) -> FieldPoint<FE> {
    //     FieldPoint { x: (), y: (), a: (), b: () }
    // }

    // fn scalar_mul(scalar: i32, point: FieldPoint<FE>) -> FieldPoint<FE> {
    //     FieldPoint { x: (), y: (), a: (), b: () }
    // }
}





#[test]
fn test_create_ec_point() {
    let new_field = Field::new(11);

    let a = FieldElement::new(5, new_field.clone());

    let b = FieldElement::new(12, new_field.clone());

    let new_curve = <EC<Field, FieldElement<Field>>>::new(a, b, new_field.clone());
}