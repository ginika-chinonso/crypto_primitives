use std::ops::{Add, BitXor, Mul};

#[derive(Debug, Clone, Copy)]
pub struct Field {
    modulus: usize,
}

#[derive(Debug, Copy, Clone)]
pub struct FieldElement {
    value: usize,
    field: Field,
}

impl Field {
    fn new(modulus: usize) -> Field {
        Field { modulus }
    }

    fn field_element(&self, value: usize) -> FieldElement {
        FieldElement::new(value, *self)
    }

    fn modulus(&self) -> usize {
        self.modulus
    }

    fn additive_identity(&self) -> FieldElement {
        FieldElement {
            value: 0,
            field: *self,
        }
    }

    fn multiplicative_identity(&self) -> FieldElement {
        FieldElement {
            value: 1,
            field: *self,
        }
    }
}

impl FieldElement {
    fn new(value: usize, field: Field) -> FieldElement {
        let val = value % field.clone().modulus();
        FieldElement { value: val, field }
    }

    fn field(&self) -> Field {
        self.field
    }

    fn value(&self) -> usize {
        self.value
    }

    fn additive_inv() {}

    fn multiplicative_inv() {}
}

impl Add for FieldElement {
    type Output = Self;
    fn add(self, rhs: Self) -> Self::Output {
        assert!(
            self.field.modulus() == rhs.field.modulus(),
            "Elements not in the same field"
        );
        let val = (self.value + rhs.value) % self.field.modulus();
        FieldElement {
            value: val,
            field: self.field,
        }
    }
}

impl Mul for FieldElement {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        assert!(self.field == rhs.field, "Elements must be in same field");
        FieldElement::new((self.value * rhs.value) % self.field.modulus(), self.field)
    }
}

// impl BitXor for FieldElement {
//     type Output = usize;

//     fn bitxor(self, rhs: usize) -> Self::Output {
//         todo!()
//     }
// }

impl PartialEq for FieldElement {
    fn eq(&self, other: &Self) -> bool {
        self.value == other.value && self.field == other.field
    }
}

impl PartialEq for Field {
    fn eq(&self, other: &Self) -> bool {
        self.modulus == other.modulus
    }
}
