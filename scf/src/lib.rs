pub mod molecular_orbitals;
pub mod utils;

use crate::{molecular_orbitals::MolecularOrbitals, utils::diis};
use basis::contracted_gaussian::ContractedGaussian;
use basis::electron_tensor::ElectronTensor;
use basis::BasisFunction;
use basis_set::atom::Atom;
use log::{debug, info, trace, warn};
use nalgebra::{DMatrix, DVector};
use std::{collections::VecDeque, time::Instant};

/// A struct containing the results of a Hartree-Fock calculation.
///
/// # Fields
///
/// * `orbitals`: `MolecularOrbitals` - Molecular orbitals obtained from the calculation.
/// * `orbital_energies`: `DVector<f64>` - Eigenvalues of the Fock matrix.
/// * `density_matrix`: `DMatrix<f64>` - Electron density matrix obtained from the calculation.
/// * `coefficient_matrix`: `DMatrix<f64>` - Coefficient matrix obtained from the calculation.
/// * `electronic_energy`: `f64` - Electronic energy obtained from the calculation.
/// * `total_energy`: `f64` - Total energy obtained from the calculation.
/// * `iterations`: `usize` - Number of iterations used to obtain the result.
/// * `n_electrons`: `usize` - Number of electrons in the system.
/// * `n_basis`: `usize` - Number of basis functions used in the calculation.
pub struct HartreeFockResult {
    pub orbitals: MolecularOrbitals,
    pub orbital_energies: DVector<f64>,
    pub density_matrix: DMatrix<f64>,
    pub coefficient_matrix: DMatrix<f64>,
    pub electronic_energy: f64,
    pub total_energy: f64,
    pub iterations: usize,
    pub n_electrons: usize,
    pub n_basis: usize,
}

/// An enum representing possible errors that may arise during a Hartree-Fock calculation.
#[derive(Copy, Clone, Debug)]
pub enum HartreeFockError {
    /// The Hartree-Fock calculation did not converge.
    NotConverged,
    /// The Hartree-Fock calculation diverged.
    Diverged,
}

pub trait SelfConsistentField {
    /// Perform a restricted Hartree-Fock calculation to approximate the electronic structure of a molecule.
    ///
    /// # Arguments
    ///
    /// * `max_iters` - Maximum number of iterations to be used in the calculation.
    /// * `epsilon` - Desired accuracy of the calculation.
    /// * `molecule_charge` - Total charge of the molecule.
    ///
    /// # Returns
    ///
    /// Returns an `Option<HartreeFockResult>`, where `HartreeFockResult` is a struct containing the HF solution.
    /// If the calculation fails to converge within the specified maximum number of iterations, `None` is returned.
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
        let n_electrons = (atoms
            .iter()
            .fold(0, |acc, atom| acc + atom.electron_count()) as i32
            - molecule_charge) as usize;

        info!("Starting SCF iteration with {n_basis} basis functions");

        if n_electrons % 2 != 0 {
            warn!("restricted hartree fock only works properly if all orbitals are fully occupied, but the specified molecule has an uneven amount of electrons");
        }

        let nuclear_repulsion = {
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
        trace!("overlap: {overlap:0.3}");
        let kinetic = utils::hermitian(n_basis, |i, j| {
            ContractedGaussian::kinetic_int(&basis[i], &basis[j])
        });
        trace!("kinetic: {kinetic:0.3}");
        let nuclear = utils::hermitian(n_basis, |i, j| {
            ContractedGaussian::nuclear_attraction_int(&basis[i], &basis[j], &point_charges)
        });
        trace!("nuclear attraction: {nuclear:0.3}");

        let start = Instant::now();
        let multi = ElectronTensor::from_basis(&basis);
        info!(
            "Multi-Electron tensor formation took {:.1?}",
            start.elapsed()
        );

        let core_hamiltonian = kinetic + nuclear;
        let transform = {
            let (u, _) = utils::eigs(overlap.clone());
            let diagonal_matrix = &u.transpose() * (&overlap * &u);

            let diagonal_inv_sqrt =
                DMatrix::from_diagonal(&diagonal_matrix.map_diagonal(|f| f.sqrt().recip()));
            &u * (diagonal_inv_sqrt * &u.transpose())
        };

        // TODO: check if this is correct, as energies in the first
        //  iteration of SCF-Cycle is a bit far off sometimes
        let mut density = {
            let hamiltonian_eht = utils::hermitian(n_basis, |i, j| {
                0.875 * overlap[(i, j)] * (core_hamiltonian[(i, i)] + core_hamiltonian[(j, j)])
            });

            let transformed = &transform.transpose() * (&hamiltonian_eht * &transform);
            let (coeffs_prime, _orbital_energies) = utils::sorted_eigs(transformed);
            let coeffs = &transform * coeffs_prime;

            utils::hermitian(n_basis, |i, j| {
                2.0 * (0..n_electrons / 2).fold(0.0, |acc, k| acc + coeffs[(i, k)] * coeffs[(j, k)])
            })
        };

        let mut previous_focks = VecDeque::new();
        let mut previous_erros = VecDeque::new();

        // precompute multi[(i, j, x, y)] - 0.5 * multi[(i, x, j, y)]
        let mut electron_terms = vec![0.0; n_basis.pow(4)];
        for j in 0..n_basis {
            for i in 0..n_basis {
                for x in 0..n_basis {
                    for y in 0..n_basis {
                        electron_terms[j * n_basis.pow(3) + i * n_basis.pow(2) + y * n_basis + x] =
                            multi[(i, j, x, y)] - 0.5 * multi[(i, x, j, y)];
                    }
                }
            }
        }

        let start = Instant::now();
        for iter in 0..=max_iters {
            let guess = utils::hermitian(n_basis, |i, j| {
                (0..n_basis).fold(0.0, |acc, y| {
                    acc + (0..n_basis).fold(0.0, |acc, x| {
                        acc + density[(x, y)]
                            * electron_terms
                                [j * n_basis.pow(3) + i * n_basis.pow(2) + y * n_basis + x]
                    })
                })
            });

            let fock = &core_hamiltonian + &guess;

            // DIIS
            // e_i = FDS - SDF
            let diis_error_estimate = &fock * &density * &overlap - &overlap * &density * &fock;
            let error = diis_error_estimate.abs().max();

            previous_erros.push_back(diis_error_estimate);
            previous_focks.push_back(fock);

            if previous_erros.len() > 16 {
                let _ = previous_erros.pop_front();
                let _ = previous_focks.pop_front();
            }

            let fock = if previous_focks.len() < 5 {
                previous_focks.back().unwrap().clone()
            } else {
                diis(&previous_erros, &previous_focks)
                    .unwrap_or_else(|| previous_focks.back().unwrap().clone())
            };
            // DIIS end

            let fock_prime = &transform.transpose() * (&fock * &transform);
            let (coeffs_prime, orbital_energies) = utils::sorted_eigs(fock_prime);
            let coeffs = &transform * &coeffs_prime;

            let new_density = utils::hermitian(n_basis, |i, j| {
                2.0 * (0..n_electrons / 2).fold(0.0, |acc, k| acc + coeffs[(i, k)] * coeffs[(j, k)])
            });

            const F: f64 = 1.0;
            let new_density = F * &new_density + (1.0 - F) * &density;

            density = new_density;

            let electronic_energy = 0.5 * (&density * (2.0 * &core_hamiltonian + &guess)).trace();
            if error < epsilon || iter == max_iters {
                let hf_energy = electronic_energy + nuclear_repulsion;
                let energies = orbital_energies.into_iter().collect::<Vec<_>>();
                let pad = (0..35).map(|_| '-').collect::<String>();

                info!("+ {pad} SCF-Routine Finished {pad} +",);
                info!("SCF took {:0.4?} to converge ({iter})", start.elapsed());
                info!("Electronic Energy: {electronic_energy:0.3}");
                info!("Nuclear Repulsion Energy: {nuclear_repulsion:0.3}");
                info!("Hartree-Fock Energy: {hf_energy:0.3}",);
                info!("Orbital Energies: {energies:0.3?}");
                info!("+ {pad} SCF-Routine Finished {pad} +",);

                return Some(HartreeFockResult {
                    orbitals: MolecularOrbitals::new(basis, &coeffs),
                    orbital_energies,
                    density_matrix: density,
                    coefficient_matrix: coeffs,
                    electronic_energy,
                    total_energy: nuclear_repulsion + electronic_energy,
                    iterations: iter,
                    n_electrons,
                    n_basis,
                });
            } else if !error.is_normal() {
                return None;
            } else {
                debug!(
                    "Iteration {iter}: max error (loss): {error:0.5e}, energy: {:0.5}, Interp: {F:0.3}",
                    electronic_energy + nuclear_repulsion
                );
            }
        }
        None
    }
}

// STO-3G with dist = 1.6a_0: https://i.imgur.com/9lXD7N7.png
// 6-31G  with dist = 1.6a_0: https://i.imgur.com/b2wa0Gd.png
#[cfg(test)]
mod tests {
    use super::*;
    use crate::utils::hermitian;
    use basis::contracted_gaussian::ContractedGaussian;
    use basis::primitives::GaussianPrimitive;
    use basis::PointCharge;
    use basis_set::periodic_table::AtomType;
    use nalgebra::base::Vector3;

    #[test]
    pub fn test_integrals() {
        let basis = &basis_set::basis_sets::BASIS_STO_3G;
        let molecule = [
            basis.get(Vector3::zeros(), AtomType::Hydrogen),
            basis.get(Vector3::new(0.0, 0.0, 1.2), AtomType::Hydrogen),
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

        let prim = GaussianPrimitive::new([0; 3], 1.0, 1.0);
        println!(
            "{:?}",
            prim.coefficient.powi(2)
                * GaussianPrimitive::_nuclear_attraction(
                    &prim,
                    &prim,
                    &Vector3::zeros(),
                    &Vector3::zeros(),
                    &PointCharge {
                        position: Vector3::new(0.0, 1.0, 0.0),
                        charge: 1.0
                    }
                )
        );

        print!("Overlap: {}", overlap);
        print!("Kinetic: {}", kinetic);
        print!("Nuclear: {}", nuclear);
    }
}
