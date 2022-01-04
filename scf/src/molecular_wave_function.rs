use basis::contracted_gaussian::ContractedGaussian;
use basis::BasisFunction;
use nalgebra::{DMatrix, Vector3};

#[derive(Debug)]
pub struct MolecularWaveFunction {
    basis_functions: Vec<ContractedGaussian>,
    coeff_matrix: DMatrix<f64>,
}

impl MolecularWaveFunction {
    const MIN_COEFFICIENT_MAGNITUDE: f64 = 1e-6;

    pub fn new(basis_functions: Vec<ContractedGaussian>, coeff_matrix: DMatrix<f64>) -> Self {
        Self {
            basis_functions,
            coeff_matrix,
        }
    }

    pub fn evaluate(&self, at: Vector3<f64>, energy_level: usize) -> f64 {
        (0..self.basis_functions.len()).fold(0.0, |acc, i| {
            let coeff = self.coeff_matrix[(i, energy_level)];
            const EPS: f64 = 5e-2;

            acc + if coeff.abs() < EPS {
                0.0
            } else {
                coeff * self.basis_functions[i].evaluate(&at)
            }
        })
    }

    pub fn basis_functions(&self) -> &Vec<ContractedGaussian> {
        &self.basis_functions
    }
}
