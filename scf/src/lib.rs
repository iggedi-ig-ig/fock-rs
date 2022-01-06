pub mod electron_tensor;
pub mod molecular_wave_function;
pub mod utils;

use crate::electron_tensor::ElectronTensor;
use crate::molecular_wave_function::MolecularWaveFunction;
use basis::contracted_gaussian::ContractedGaussian;
use basis::BasisFunction;
use basis_set::atom::Atom;
use nalgebra::{DMatrix, DVector};
use std::time::Instant;

#[derive(Debug)]
pub struct HartreeFockResult {
    pub wave_function: MolecularWaveFunction,
    pub orbital_energies: DVector<f64>,
    pub density_matrix: DMatrix<f64>,
    pub electronic_energy: f64,
    pub total_energy: f64,
    pub iterations: usize,
    pub n_electrons: usize,
}

#[derive(Copy, Clone, Debug)]
pub enum HartreeFockError {
    NotConverged,
    Diverged,
}

pub trait SelfConsistentField {
    fn try_scf(
        &self,
        max_iters: usize,
        epsilon: f64,
        molecule_charge: i32,
    ) -> Option<HartreeFockResult>;
}

impl<T> SelfConsistentField for T
where
    T: IntoIterator<Item = Atom> + Clone,
{
    fn try_scf(
        &self,
        max_iters: usize,
        epsilon: f64,
        molecule_charge: i32,
    ) -> Option<HartreeFockResult> {
        let atoms = self.clone().into_iter().collect::<Vec<Atom>>();
        let basis = atoms
            .iter()
            .flat_map(|atom| atom.basis().iter().cloned())
            .collect::<Vec<_>>();
        let point_charges = atoms
            .iter()
            .map(|atom| atom.point_charge())
            .collect::<Vec<_>>();

        let n_atoms = atoms.len();
        let n_basis = basis.len();
        let n_elecs = (atoms
            .iter()
            .fold(0, |acc, atom| acc + atom.valence_electrons()) as i32
            - molecule_charge) as usize;

        let nuclear_repulsion_energy = {
            let point_charges = &point_charges;

            (0..n_atoms)
                .flat_map(|i| {
                    (i + 1..n_atoms).map(move |j| {
                        let atom_a = &point_charges[i];
                        let atom_b = &point_charges[j];

                        let diff = atom_b.position - atom_a.position;
                        atom_a.charge * atom_b.charge / diff.norm()
                    })
                })
                .sum::<f64>()
        };

        let overlap = utils::hermitian(n_basis, |i, j| {
            ContractedGaussian::overlap_int(&basis[i], &basis[j])
        });
        let kinetic = utils::hermitian(n_basis, |i, j| {
            ContractedGaussian::kinetic_int(&basis[i], &basis[j])
        });
        let nuclear = utils::hermitian(n_basis, |i, j| {
            ContractedGaussian::nuclear_attraction_int(&basis[i], &basis[j], &point_charges)
        });

        let start = Instant::now();
        let multi = ElectronTensor::from_basis(&basis);
        println!(
            "\rMulti-Electron tensor formation took {:0.4?}",
            start.elapsed()
        );

        let core_hamiltonian = kinetic + nuclear;
        let transformation = {
            let (unitary, _) = utils::sorted_eigs(overlap.clone());
            let diagonalized = &unitary.transpose() * (&overlap * &unitary);
            let diagonal = DMatrix::from_diagonal(&diagonalized.map_diagonal(|f| f.sqrt().recip()));
            &unitary * (&diagonal * &unitary.transpose())
        };
        let transform = |m: &DMatrix<f64>| &transformation.transpose() * (m * &transformation);
        let transform_inv = |m: &DMatrix<f64>| &transformation * m;

        let mut density = DMatrix::from_element(n_basis, n_basis, 0.0f64);
        for iteration in 0..max_iters {
            let guess = {
                let multi = &multi;
                let density = &density;

                DMatrix::from_fn(n_basis, n_basis, |i, j| {
                    (0..n_basis)
                        .flat_map(|x| {
                            (0..n_basis).map(move |y| {
                                density[(x, y)] * (multi[(i, j, x, y)] - 0.5 * multi[(i, x, y, j)])
                            })
                        })
                        .sum::<f64>()
                })
            };

            let fock = &core_hamiltonian + &guess;
            let fock_prime = transform(&fock);

            let (coeffs_prime, orbital_energies) = utils::sorted_eigs(fock_prime);
            let coeffs = transform_inv(&coeffs_prime);

            let new_density = DMatrix::from_fn(n_basis, n_basis, |i, j| {
                2.0 * (0..(n_elecs + 1) / 2)
                    .fold(0.0, |acc, k| acc + coeffs[(i, k)] * coeffs[(j, k)])
            });

            let density_rms = density
                .zip_fold(&new_density, 0.0, |acc, new, old| acc + (new - old).powi(2))
                .sqrt()
                / n_basis as f64;

            if density_rms < epsilon {
                let electronic_energy = 0.5 * (&new_density * (&core_hamiltonian + &fock)).trace();

                return Some(HartreeFockResult {
                    wave_function: MolecularWaveFunction::new(basis, coeffs),
                    orbital_energies,
                    density_matrix: new_density,
                    electronic_energy,
                    total_energy: nuclear_repulsion_energy + electronic_energy,
                    iterations: iteration,
                    n_electrons: n_elecs,
                });
            } else if !density_rms.is_normal() {
                return None;
            }

            density = new_density;
        }
        None
    }
}

// STO-3G with dist = 1.6a_0: https://i.imgur.com/9lXD7N7.png
// 6-31G  with dist = 1.6a_0: https://i.imgur.com/b2wa0Gd.png
#[test]
pub fn test_integrals() {
    use crate::utils::hermitian;
    use basis_set::periodic_table::AtomType;
    use nalgebra::base::Vector3;

    let basis = &basis_set::basis_sets::BASIS_STO_3G;
    let molecule = [
        basis.get(Vector3::zeros(), AtomType::Hydrogen),
        basis.get(Vector3::new(0.0, 0.0, 1.6), AtomType::Hydrogen),
    ];

    let basis = molecule
        .iter()
        .flat_map(|atom| atom.basis().clone())
        .collect::<Vec<_>>();
    let point_charges = molecule
        .iter()
        .map(|atom| atom.point_charge())
        .collect::<Vec<_>>();
    let n = basis.len();

    let overlap = hermitian(n, |i, j| {
        ContractedGaussian::overlap_int(&basis[i], &basis[j])
    });
    let kinetic = hermitian(n, |i, j| {
        ContractedGaussian::kinetic_int(&basis[i], &basis[j])
    });
    let nuclear = hermitian(n, |i, j| {
        ContractedGaussian::nuclear_attraction_int(&basis[i], &basis[j], &point_charges)
    });

    print!("Overlap: {:0.5}", overlap);
    print!("Kinetic: {:0.5}", kinetic);
    print!("Nuclear: {:0.5}", nuclear);
}
