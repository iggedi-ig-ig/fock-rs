#![allow(dead_code)]

use basis_set::atom::Atom;
use basis_set::periodic_table::AtomType;
use basis_set::periodic_table::AtomType::*;
use basis_set::BasisSet;
use nalgebra::Vector3;

macro_rules! bohr {
    ($num:expr) => {
        $num * 1.89
    };
}

macro_rules! define_molecules {
    ($($molecule:ident = [$($atom:path = ($x:expr, $y:expr, $z:expr)),*]),*) => {
        $(
        pub const $molecule: &crate::molecules::MoleculeBlueprint = &crate::molecules::MoleculeBlueprint {
            atoms: &[
                $(
                    crate::molecules::AtomBlueprint::new(nalgebra::Vector3::new($x, $y, $z), $atom),
                )*
            ]
        };
        )*
    }
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

define_molecules! {
    HYDROGEN = [Hydrogen = (-0.8, 0.0, 0.0), Hydrogen = (0.8, 0.0, 0.0)],
    HELIUM_HYDRIDE = [Hydrogen = (-0.7, 0.0, 0.0), Helium = (0.7, 0.0, 0.0)],
    WATER = [
        Hydrogen = (-0.78 * bohr!(0.96), -0.62 * bohr!(0.96), 0.0),
        Oxygen = (0.0, 0.0, 0.0),
        Hydrogen = (0.78 * bohr!(0.96), -0.62 * bohr!(0.96), 0.0)
    ],
    NITRITE = [
        Oxygen = (-0.881 * bohr!(1.24), -0.4732 * bohr!(1.24), 0.0),
        Nitrogen = (0.0, 0.0, 0.0),
        Oxygen = (0.881 * bohr!(1.24), -0.4732 * bohr!(1.24), 0.0)
    ],
    ETHENE = [
        Hydrogen = (-1.25 - 0.866 * 2.0, 0.0, -1.0),
        Hydrogen = (-1.25 - 0.866 * 2.0, 0.0, 1.0),
        Carbon = (-1.25, 0.0, 0.0),
        Carbon = (1.25, 0.0, 0.0),
        Hydrogen = (1.25 + 0.866 * 2.0, 0.0, -1.0),
        Hydrogen = (1.25 + 0.866 * 2.0, 0.0, 1.0)
    ],
    ETHYLENE = [
        Hydrogen = (-3.125, 0.0, 0.0),
        Carbon = (-1.125, 0.0, 0.0),
        Carbon = (1.125, 0.0, 0.0),
        Hydrogen = (3.125, 0.0, 0.0)
    ],
    BENZENE = [
        Carbon = (0.00000e0, 0.0, 2.62672e0),
        Hydrogen = (0.00000e0, 0.0, 4.68652e0),
        Carbon = (2.27481e0, 0.0, 1.31336e0),
        Hydrogen = (4.05865e0, 0.0, 2.34326e0),
        Carbon = (2.27481e0, 0.0, -1.31336e0),
        Hydrogen = (4.05865e0, 0.0, -2.34326e0),
        Carbon = (-8.44817e-16, 0.0, -2.62672e0),
        Hydrogen = (-1.50730e-15, 0.0, -4.68652e0),
        Carbon = (-2.27481e0, 0.0, -1.31336e0),
        Hydrogen = (-4.05865e0, 0.0, -2.34326e0),
        Carbon = (-2.27481e0, 0.0, 1.31336e0),
        Hydrogen = (-4.05865e0, 0.0, 2.34326e0)
    ],
    METHANE = [
        Carbon = (0.0, 0.0, 0.0),
        Hydrogen = (0.0, bohr!(1.09), 0.0),
        Hydrogen = (0.0, -0.3338 * bohr!(1.09), bohr!(1.09)),
        Hydrogen = (0.866 * bohr!(1.09), -0.3338 * bohr!(1.09), -0.5 * bohr!(1.09))
    ],
    AMMONIA = [
        Nitrogen = (0.0, 0.0, 0.0),
        Hydrogen = (0.0, -0.3338 * bohr!(1.017), bohr!(1.017)),
        Hydrogen = (0.866 * bohr!(1.02), -0.3338 * bohr!(1.02), -0.5 * bohr!(1.02)),
        Hydrogen = (-0.866 * bohr!(1.02), -0.3338 * bohr!(1.02), -0.5 * bohr!(1.02))
    ],
    NAPHTALENE = [
        Carbon = (2.627, 0.0, 2.627),
        Hydrogen = (4.687, 0.0, 4.687),
        Carbon = (4.902, 0.0, 1.314),
        Hydrogen = (8.746, 0.0, 2.344),
        Carbon = (4.902, 0.0, -1.314),
        Hydrogen = (8.746, 0.0, -2.344),
        Carbon = (2.627, 0.0, -2.627),
        Hydrogen = (4.687, 0.0, -4.687),
        Carbon = (-2.627, 0.0, 2.627),
        Hydrogen = (-4.687, 0.0, 4.687),
        Carbon = (-4.902, 0.0, 1.314),
        Hydrogen = (-8.746, 0.0, 2.344),
        Carbon = (-4.902, 0.0, -1.314),
        Hydrogen = (-8.746, 0.0, -2.344),
        Carbon = (-2.627, 0.0, -2.627),
        Hydrogen = (-4.687, 0.0, -4.687),
        Carbon = (0.0, 0.0, 1.314),
        Carbon = (0.0, 0.0, -1.314)
    ]
}
