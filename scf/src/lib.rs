pub mod electron_tensor;
pub mod molecular_wave_function;
pub mod utils;

use crate::electron_tensor::ElectronRepulsionTensor;
use crate::molecular_wave_function::MolecularWaveFunction;
use crate::utils::hermitian;
use basis::contracted_gaussian::ContractedGaussian;
use basis::BasisFunction;
use basis_set::atom::Atom;
use basis_set::periodic_table::AtomType;
use nalgebra::{DMatrix, DVector, Vector3};

#[derive(Debug)]
pub struct HartreeFockResult {
    pub wave_function: MolecularWaveFunction,
    pub orbital_energies: DVector<f64>,
    pub density_matrix: DMatrix<f64>,
    pub electronic_energy: f64,
    pub total_energy: f64,
    pub iterations: usize,
}

#[derive(Copy, Clone, Debug)]
pub enum HartreeFockError {
    NotConverged,
    Diverged,
}

pub trait SelfConsistentField {
    fn try_scf(&self, n_elecs: usize, max_iters: usize, epsilon: f64) -> Option<HartreeFockResult>;
}

impl<T> SelfConsistentField for T
where
    T: IntoIterator<Item = Atom> + Clone,
{
    fn try_scf(&self, n_elecs: usize, max_iters: usize, epsilon: f64) -> Option<HartreeFockResult> {
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

        // println!("Num atoms: {}", n_atoms);
        // println!("Num electrons: {}", n_elecs);
        // println!("Num basis functions: {}", n_basis);

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
        let multi = ElectronRepulsionTensor::from_fn(n_basis, |i, j, k, l| {
            ContractedGaussian::electron_repulsion_int(&basis[i], &basis[j], &basis[k], &basis[l])
        });

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
                                density[(x, y)] * (multi[(i, j, x, y)] - 0.5 * multi[(i, y, x, j)])
                            })
                        })
                        .sum::<f64>()
                })
            };

            let fock = &core_hamiltonian + &guess;
            let fock_prime = transform(&fock);

            let (coeffs_prime, orbital_energies) = utils::sorted_eigs(fock_prime);
            let coeffs = &transformation * coeffs_prime;

            let new_density = DMatrix::from_fn(n_basis, n_basis, |i, j| {
                2.0 * (0..n_elecs / 2).fold(0.0, |acc, k| acc + coeffs[(i, k)] * coeffs[(j, k)])
            });

            let density_rms = density
                .zip_fold(&new_density, 0.0, |acc, new, old| acc + (new - old).powi(2))
                .sqrt()
                / n_basis as f64;

            // println!("Iteration {}: rms: {:0.5e}", iteration, density_rms);
            if density_rms < epsilon {
                // println!("Converged!");
                let electronic_energy = 0.5 * (&new_density * (&core_hamiltonian + &fock)).trace();

                return Some(HartreeFockResult {
                    wave_function: MolecularWaveFunction::new(basis, coeffs),
                    orbital_energies,
                    density_matrix: new_density,
                    electronic_energy,
                    total_energy: nuclear_repulsion_energy + electronic_energy,
                    iterations: iteration,
                });
            }

            density = new_density;
        }
        None
    }
}

#[test]
pub fn test_plot_dist() {
    fn test(dist: f64) {
        let basis = &basis_set::basis_sets::BASIS_6_31G;
        let molecule = [
            basis
                .get(Vector3::new(0.0, 0.0, -dist * 0.5), AtomType::Helium, 0)
                .unwrap(),
            basis
                .get(Vector3::new(0.0, 0.0, dist * 0.5), AtomType::Hydrogen, -1)
                .unwrap(),
        ];

        let result = molecule.try_scf(2, 1000, 1e-6).unwrap();
        print!("{:0.4}, ", result.total_energy);
    }

    print!("data = [");
    for dist in (1..250).map(|i| i as f64 / 250.0 * 2.5) {
        test(dist);
    }
    print!("]")
}

#[test]
pub fn test_helium_hydride_cation() {
    const DIST: f64 = 1.4;

    let basis = &basis_set::basis_sets::BASIS_6_31G;
    let molecule = [
        basis
            .get(Vector3::new(0.0, 0.0, -DIST * 0.5), AtomType::Helium, 0)
            .unwrap(),
        basis
            .get(Vector3::new(0.0, 0.0, DIST * 0.5), AtomType::Hydrogen, -1)
            .unwrap(),
    ];

    if let Some(result) = molecule.try_scf(2, 1000, 1e-5) {
        println!("SCF-Cycle converged in {} iterations", result.iterations);
        print!(
            "Orbital (coefficient) matrix: {:0.5}",
            result.wave_function.coeff_matrix()
        );
        print!("Density matrix: {:0.5}", result.density_matrix);
        print!(
            "Orbital energies: {:0.5}",
            result.orbital_energies.transpose()
        );
        println!("Total energy: {:0.4}", result.total_energy);
        println!("Electronic energy: {:0.4}", result.electronic_energy);
    }
}

// STO-3G with dist = 1.6a_0: https://i.imgur.com/9lXD7N7.png
// 6-31G  with dist = 1.6a_0: https://i.imgur.com/b2wa0Gd.png
#[test]
pub fn test_integrals() {
    let basis = &basis_set::basis_sets::BASIS_STO_3G;
    let molecule = [
        basis.get(Vector3::zeros(), AtomType::Hydrogen, 0).unwrap(),
        basis
            .get(Vector3::new(0.0, 0.0, 1.6), AtomType::Hydrogen, 0)
            .unwrap(),
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
