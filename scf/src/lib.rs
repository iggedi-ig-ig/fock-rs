pub mod electron_tensor;
pub mod molecular_orbitals;
pub mod utils;

use crate::electron_tensor::ElectronTensor;
use crate::molecular_orbitals::MolecularOrbitals;
use crate::utils::hermitian;
use basis::contracted_gaussian::ContractedGaussian;
use basis::BasisFunction;
use basis_set::atom::Atom;
use log::{debug, info, warn};
use nalgebra::{DMatrix, DVector};
use std::time::Instant;

#[derive(Debug)]
pub struct HartreeFockResult {
    pub orbitals: MolecularOrbitals,
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
            .fold(0, |acc, atom| acc + atom.electron_count()) as i32
            - molecule_charge) as usize;

        if n_elecs % 2 != 0 {
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
        debug!("overlap: {overlap:0.5}");
        let kinetic = utils::hermitian(n_basis, |i, j| {
            ContractedGaussian::kinetic_int(&basis[i], &basis[j])
        });
        debug!("kinetic: {kinetic:0.5}");
        let nuclear = utils::hermitian(n_basis, |i, j| {
            ContractedGaussian::nuclear_attraction_int(&basis[i], &basis[j], &point_charges)
        });
        debug!("nuclear attraction: {nuclear:0.5}");

        let start = Instant::now();
        let multi = ElectronTensor::from_basis(&basis);
        info!(
            "Multi-Electron tensor formation took {:0.4?}",
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

        let start = Instant::now();
        let mut density = DMatrix::from_element(n_basis, n_basis, 0.0f64);
        for iter in 0..max_iters {
            let guess = hermitian(n_basis, |i, j| {
                (0..n_basis).fold(0.0, |acc, x| {
                    acc + (0..n_basis).fold(0.0, |acc, y| {
                        acc + density[(x, y)] * (multi[(i, j, y, x)] - 0.5 * multi[(i, x, y, j)])
                    })
                })
            });

            let fock = &core_hamiltonian + &guess;
            let fock_prime = &transform.transpose() * (&fock * &transform);

            let (coeffs_prime, orbital_energies) = utils::sorted_eigs(fock_prime);
            let coeffs = &transform * &coeffs_prime;

            let new_density = utils::hermitian(n_basis, |i, j| {
                2.0 * (0..n_elecs / 2).fold(0.0, |acc, k| acc + coeffs[(i, k)] * coeffs[(j, k)])
            });

            let density_rms = density
                .zip_fold(&new_density, 0.0, |acc, new, old| acc + (new - old).powi(2))
                .sqrt()
                / n_basis as f64;

            let electronic_energy = 0.5 * (&density * (2.0 * &core_hamiltonian + &guess)).trace();
            if density_rms < epsilon {
                let pad = (0..35).map(|_| '-').collect::<String>();
                println!("+ {pad} SCF-Routine Finished {pad} +",);
                println!(
                    "SCF took {:0.4?} to converge ({iter} iters)",
                    start.elapsed(),
                );

                println!("Electronic Energy: {electronic_energy:0.4}");
                println!("Nuclear Repulsion Energy: {nuclear_repulsion:0.4}");
                println!(
                    "Total Hartree-Fock Energy: {:0.4}",
                    electronic_energy + nuclear_repulsion
                );
                println!("Final Density Matrix: {new_density:0.5}");
                println!("Final Coefficient Matrix: {coeffs:0.5}");
                println!(
                    "Final Orbital Energies: {:0.5}",
                    orbital_energies.transpose()
                );
                println!("+ {pad} SCF-Routine Finished {pad} +",);

                return Some(HartreeFockResult {
                    orbitals: MolecularOrbitals::new(basis, coeffs),
                    orbital_energies,
                    density_matrix: new_density,
                    electronic_energy,
                    total_energy: nuclear_repulsion + electronic_energy,
                    iterations: iter,
                    n_electrons: n_elecs,
                });
            } else if !density_rms.is_normal() {
                return None;
            } else {
                debug!(
                    "Iteration {iter}: density rms: {density_rms:0.5e}, energy: {:0.5}",
                    electronic_energy + nuclear_repulsion
                );
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
