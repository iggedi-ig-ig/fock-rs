use crate::AtomType;
use basis::contracted_gaussian::ContractedGaussian;
use basis::PointCharge;
use nalgebra::Vector3;
use std::fmt::{Display, Formatter};

#[derive(Clone, Debug)]
pub struct Atom {
    position: Vector3<f64>,
    atom_type: AtomType,
    /// this is the amount of exta electrons this atom has
    ///
    /// example:
    ///     - O^2- has an electron balance of -2
    ///     - H+   has an electron balance of +1
    /// etc.
    electron_balance: i32,
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
            electron_balance: 0,
            basis,
        }
    }

    pub fn new_ion(
        position: Vector3<f64>,
        atom_type: AtomType,
        basis: Vec<ContractedGaussian>,
        electron_balance: i32,
    ) -> Self {
        Self {
            position,
            atom_type,
            electron_balance,
            basis,
        }
    }

    pub fn position(&self) -> Vector3<f64> {
        self.position
    }
    pub fn basis(&self) -> &Vec<ContractedGaussian> {
        &self.basis
    }
    pub fn num_electrons(&self) -> usize {
        ((self.atom_type as usize) as i32 + self.electron_balance) as usize
    }
    pub fn atom_type(&self) -> AtomType {
        self.atom_type
    }
    pub fn electron_balance(&self) -> i32 {
        self.electron_balance
    }
    pub fn point_charge(&self) -> PointCharge {
        PointCharge {
            position: self.position(),
            charge: (self.atom_type as usize) as f64,
        }
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
