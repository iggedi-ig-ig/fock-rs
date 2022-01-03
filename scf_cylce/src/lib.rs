pub mod molecular_wave_function;

use crate::molecular_wave_function::MolecularWaveFunction;
use basis::BasisFunction;
use nalgebra::{DMatrix, DVector};

pub struct HartreeFockResult {
    pub wave_function: MolecularWaveFunction,
    pub orbital_energies: DVector<f64>,
    pub density_matrix: DMatrix<f64>,
    pub electronic_energy: f64,
    pub total_energy: f64,
}

#[derive(Copy, Clone, Debug)]
pub enum HartreeFockError {
    NotConverged,
    Diverged,
}

pub trait SelfConsistentField {
    fn try_scf(&self, max_iters: usize, epsilon: f64) -> Option<HartreeFockResult>;
}
