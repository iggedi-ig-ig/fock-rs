#![allow(dead_code)]

use basis_set::atom::Atom;
use basis_set::periodic_table::AtomType;
use basis_set::periodic_table::AtomType::{Carbon, Hydrogen, Nitrogen, Oxygen};
use basis_set::BasisSet;
use nalgebra::Vector3;

#[macro_export]
macro_rules! bohr {
    ($num:expr) => {
        $num * 1.89
    };
}

#[derive(Copy, Clone)]
pub struct AtomBlueprint {
    position: Vector3<f64>,
    atom_type: AtomType,
}

impl AtomBlueprint {
    pub const fn new(position: Vector3<f64>, atom_type: AtomType) -> Self {
        Self {
            position,
            atom_type,
        }
    }

    pub fn build(self, basis_set: &BasisSet) -> Atom {
        basis_set.get(self.position, self.atom_type)
    }
}

#[derive(Default)]
pub struct MoleculeBlueprint<'a> {
    atoms: &'a [AtomBlueprint],
}

impl<'a> MoleculeBlueprint<'a> {
    pub fn build(&self, basis_set: &BasisSet) -> Vec<Atom> {
        self.atoms
            .iter()
            .map(|atom| atom.build(basis_set))
            .collect()
    }
}

pub const HYDROGEN: &MoleculeBlueprint = &MoleculeBlueprint {
    atoms: &[
        AtomBlueprint::new(Vector3::new(-0.8, 0.0, 0.0), Hydrogen),
        AtomBlueprint::new(Vector3::new(0.8, 0.0, 0.0), Hydrogen),
    ],
};
pub const WATER: &MoleculeBlueprint = &MoleculeBlueprint {
    atoms: &[
        AtomBlueprint::new(
            Vector3::new(-0.78 * bohr!(0.96), -0.62 * bohr!(0.96), 0.0),
            Hydrogen,
        ),
        AtomBlueprint::new(Vector3::new(0.0, 0.0, 0.0), Oxygen),
        AtomBlueprint::new(
            Vector3::new(0.78 * bohr!(0.96), -0.62 * bohr!(0.96), 0.0),
            Hydrogen,
        ),
    ],
};
pub const NITRITE: &MoleculeBlueprint = &MoleculeBlueprint {
    atoms: &[
        AtomBlueprint::new(
            Vector3::new(-0.881 * bohr!(1.24), -0.4732 * bohr!(1.24), 0.0),
            Oxygen,
        ),
        AtomBlueprint::new(Vector3::new(0.0, 0.0, 0.0), Nitrogen),
        AtomBlueprint::new(
            Vector3::new(0.881 * bohr!(1.24), -0.4732 * bohr!(1.24), 0.0),
            Oxygen,
        ),
    ],
};
pub const ETHENE: &MoleculeBlueprint = &MoleculeBlueprint {
    atoms: &[
        AtomBlueprint::new(Vector3::new(-1.25 - 0.866 * 2.0, 0.0, -1.0), Hydrogen),
        AtomBlueprint::new(Vector3::new(-1.25 - 0.866 * 2.0, 0.0, 1.0), Hydrogen),
        AtomBlueprint::new(Vector3::new(-1.25, 0.0, 0.0), Carbon),
        AtomBlueprint::new(Vector3::new(1.25, 0.0, 0.0), Carbon),
        AtomBlueprint::new(Vector3::new(1.25 + 0.866 * 2.0, 0.0, -1.0), Hydrogen),
        AtomBlueprint::new(Vector3::new(1.25 + 0.866 * 2.0, 0.0, 1.0), Hydrogen),
    ],
};
pub const ETHYLENE: &MoleculeBlueprint = &MoleculeBlueprint {
    atoms: &[
        AtomBlueprint::new(Vector3::new(-3.125, 0.0, 0.0), Hydrogen),
        AtomBlueprint::new(Vector3::new(-1.125, 0.0, 0.0), Carbon),
        AtomBlueprint::new(Vector3::new(1.125, 0.0, 0.0), Carbon),
        AtomBlueprint::new(Vector3::new(3.125, 0.0, 0.0), Hydrogen),
    ],
};
pub const BENZENE: &MoleculeBlueprint = &MoleculeBlueprint {
    atoms: &[
        AtomBlueprint::new(Vector3::new(0.00000e0, 0.0, 2.62672e0), Carbon),
        AtomBlueprint::new(Vector3::new(0.00000e0, 0.0, 4.68652e0), Hydrogen),
        AtomBlueprint::new(Vector3::new(2.27481e0, 0.0, 1.31336e0), Carbon),
        AtomBlueprint::new(Vector3::new(4.05865e0, 0.0, 2.34326e0), Hydrogen),
        AtomBlueprint::new(Vector3::new(2.27481e0, 0.0, -1.31336e0), Carbon),
        AtomBlueprint::new(Vector3::new(4.05865e0, 0.0, -2.34326e0), Hydrogen),
        AtomBlueprint::new(Vector3::new(-8.44817e-16, 0.0, -2.62672e0), Carbon),
        AtomBlueprint::new(Vector3::new(-1.50730e-15, 0.0, -4.68652e0), Hydrogen),
        AtomBlueprint::new(Vector3::new(-2.27481e0, 0.0, -1.31336e0), Carbon),
        AtomBlueprint::new(Vector3::new(-4.05865e0, 0.0, -2.34326e0), Hydrogen),
        AtomBlueprint::new(Vector3::new(-2.27481e0, 0.0, 1.31336e0), Carbon),
        AtomBlueprint::new(Vector3::new(-4.05865e0, 0.0, 2.34326e0), Hydrogen),
    ],
};
pub const METHANE: &MoleculeBlueprint = &MoleculeBlueprint {
    atoms: &[
        AtomBlueprint::new(Vector3::new(0.0, 0.0, 0.0), Carbon),
        AtomBlueprint::new(Vector3::new(0.0, bohr!(1.09), 0.0), Hydrogen),
        AtomBlueprint::new(
            Vector3::new(0.0, -0.3338 * bohr!(1.09), bohr!(1.09)),
            Hydrogen,
        ),
        AtomBlueprint::new(
            Vector3::new(
                0.866 * bohr!(1.09),
                -0.3338 * bohr!(1.09),
                -0.5 * bohr!(1.09),
            ),
            Hydrogen,
        ),
        AtomBlueprint::new(
            Vector3::new(
                -0.866 * bohr!(1.09),
                -0.3338 * bohr!(1.09),
                -0.5 * bohr!(1.09),
            ),
            Hydrogen,
        ),
    ],
};
pub const AMMONIA: &MoleculeBlueprint = &MoleculeBlueprint {
    atoms: &[
        AtomBlueprint::new(Vector3::new(0.0, 0.0, 0.0), Nitrogen),
        AtomBlueprint::new(
            Vector3::new(0.0, -0.3338 * bohr!(1.017), bohr!(1.017)),
            Hydrogen,
        ),
        AtomBlueprint::new(
            Vector3::new(
                0.866 * bohr!(1.02),
                -0.3338 * bohr!(1.02),
                -0.5 * bohr!(1.02),
            ),
            Hydrogen,
        ),
        AtomBlueprint::new(
            Vector3::new(
                -0.866 * bohr!(1.02),
                -0.3338 * bohr!(1.02),
                -0.5 * bohr!(1.02),
            ),
            Hydrogen,
        ),
    ],
};

#[test]
pub fn gen_ammonia() {
    for angle in (0..3).map(|k| k as f32 * 2.0 * std::f32::consts::FRAC_PI_3) {
        println!("AtomBlueprint::new(Vector3::new({:0.3} * bohr!(1.02), -0.3338 * bohr!(1.02), {:0.3} * bohr!(1.02)), Hydrogen),", angle.sin(), angle.cos())
    }
}
