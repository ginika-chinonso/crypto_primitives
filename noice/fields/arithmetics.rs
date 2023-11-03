#[derive(Debug, Clone, Copy)]
pub struct Field{
    modulus: i32,
}


#[derive(Debug, Copy, Clone)]
pub struct FieldElement<F> where F: FieldTrait + Clone + Copy {
    value: i32,
    field: F
}

pub trait FieldTrait{
    fn new(modulus: i32) -> Field;
    fn modulus(self: Self) -> i32;
}


pub trait FieldElementTrait<F: FieldTrait + Clone + Copy>{
    fn new(value: i32, field: F) -> FieldElement<F>;
    fn value(self: Self) -> i32;
    fn field(self: Self) -> F;
    fn add(self: Self, rhs: FieldElement<F>) -> FieldElement<F>;
}


impl FieldTrait for Field {
    fn new(modulus: i32) -> Field {
        Field { modulus }
    }

    fn modulus(self: Self) -> i32 {
        self.modulus
    }
}


impl <F: FieldTrait + Clone + Copy> FieldElementTrait<F> for FieldElement<F> {
    fn new(value: i32, field: F) -> FieldElement<F> {
        let val = value % field.clone().modulus();
        FieldElement { value: val, field }
    }

    
    fn field(self: Self) -> F {
        self.field
    }
    
    fn value(self: Self) -> i32 {
        self.value
    }
    
    fn add(self: Self, rhs: FieldElement<F>) -> FieldElement<F> {
        assert!(self.field.modulus() == rhs.field.modulus(), "Elements not in the same field");
        let val = (self.value + rhs.value) % self.field.modulus();
        FieldElement { value: val, field: self.field }
    }
}