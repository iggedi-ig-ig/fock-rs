use crate::AtomType;
use basis::contracted_gaussian::ContractedGaussian;
use basis::PointCharge;
use nalgebra::Vector3;

#[derive(Clone, Debug)]
pub struct Atom {
    position: Vector3<f64>,
    atom_type: AtomType,
    basis: Vec<ContractedGaussian>,
}

impl Atom {
    pub fn new(
        position: Vector3<f64>,
        atom_type: AtomType,
        basis: Vec<ContractedGaussian>,
    ) -> Self {
        Self {
            position,
            atom_type,
            basis,
        }
    }

    pub fn position(&self) -> Vector3<f64> {
        self.position
    }
    pub fn basis(&self) -> &Vec<ContractedGaussian> {
        &self.basis
    }
    pub fn electron_count(&self) -> usize {
        self.atom_type as usize
    }
    pub fn atom_type(&self) -> AtomType {
        self.atom_type
    }
    pub fn point_charge(&self) -> PointCharge {
        PointCharge {
            position: self.position(),
            charge: (self.atom_type as usize) as f64,
        }
    }
}
