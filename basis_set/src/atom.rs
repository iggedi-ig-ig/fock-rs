use crate::AtomType;
use basis::contracted_gaussian::ContractedGaussian;
use nalgebra::Vector3;
use std::fmt::{Display, Formatter};

#[derive(Clone, Debug)]
pub struct Atom {
    position: Vector3<f64>,
    atom_type: AtomType,
    oxidation_state: i32,
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
            oxidation_state: 0,
            basis,
        }
    }

    pub fn new_ion(
        position: Vector3<f64>,
        atom_type: AtomType,
        basis: Vec<ContractedGaussian>,
        oxidation_state: i32,
    ) -> Self {
        Self {
            position,
            atom_type,
            oxidation_state,
            basis,
        }
    }

    pub fn position(&self) -> Vector3<f64> {
        self.position
    }
    pub fn ordinal(&self) -> usize {
        self.atom_type as usize
    }
    pub fn basis(&self) -> &[ContractedGaussian] {
        &self.basis
    }
    pub fn num_electrons(&self) -> usize {
        (self.atom_type as i32 + self.oxidation_state as i32) as usize
    }
}

impl Display for Atom {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        writeln!(
            f,
            "{:?} (electrons: {})\n\tPos: {:0.3}, {:0.3}, {:0.3}\n\tBasis Functions:",
            self.atom_type,
            self.num_electrons(),
            self.position.x,
            self.position.y,
            self.position.z
        )?;
        for (i, contracted_gaussian) in self.basis.iter().enumerate() {
            writeln!(f, "\t{}:", i)?;
            for primitive in contracted_gaussian.primitives() {
                writeln!(f, "\t\t{}", primitive)?;
            }
        }
        Ok(())
    }
}
