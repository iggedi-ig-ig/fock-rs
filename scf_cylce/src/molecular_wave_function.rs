use basis::BasisFunction;
use nalgebra::{DMatrix, Vector3};

pub struct MolecularWaveFunction<B: BasisFunction> {
    basis_functions: Vec<B>,
    coeff_matrix: DMatrix<f64>,
}

impl<B: BasisFunction> MolecularWaveFunction<B> {
    const MIN_COEFFICIENT_MAGNITUDE: f64 = 1e-6;

    pub fn evaluate(&self, at: Vector3<f64>, energy_level: usize) -> f64 {
        let coeffs = self.coeff_matrix.column(energy_level);
        coeffs
            .into_iter()
            .zip(self.basis_functions.iter())
            .filter(|(coeff, _)| coeff.abs() > Self::MIN_COEFFICIENT_MAGNITUDE)
            .map(|(coeff, basis)| coeff * basis.evaluate(&at))
            .sum::<f64>()
    }
}
