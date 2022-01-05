use basis::contracted_gaussian::ContractedGaussian;
use basis::BasisFunction;
use nalgebra::{DMatrix, Vector3};

#[derive(Debug)]
pub struct MolecularWaveFunction {
    basis_functions: Vec<ContractedGaussian>,
    coeff_matrix: DMatrix<f64>,
}

impl MolecularWaveFunction {
    const MIN_COEFFICIENT_MAGNITUDE: f64 = 0.05;

    pub fn new(basis_functions: Vec<ContractedGaussian>, coeff_matrix: DMatrix<f64>) -> Self {
        Self {
            basis_functions,
            coeff_matrix,
        }
    }

    pub fn evaluate(&self, at: Vector3<f64>, energy_level: usize) -> f64 {
        let coeffs = self.coeff_matrix.column(energy_level);
        coeffs
            .into_iter()
            .zip(self.basis_functions.iter())
            .filter(|(coeff, _)| coeff.abs() > Self::MIN_COEFFICIENT_MAGNITUDE)
            .map(|(coeff, basis)| coeff * basis.evaluate(&at))
            .sum::<f64>()
    }

    pub fn basis_functions(&self) -> &Vec<ContractedGaussian> {
        &self.basis_functions
    }
    pub fn coeff_matrix(&self) -> &DMatrix<f64> {
        &self.coeff_matrix
    }
}
