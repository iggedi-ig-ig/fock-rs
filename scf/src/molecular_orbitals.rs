use basis::contracted_gaussian::ContractedGaussian;
use basis::BasisFunction;
use nalgebra::{DMatrix, DVector, Vector3};
use std::ops::Index;

pub struct MolecularOrbitals {
    coeff_matrix: DMatrix<f64>,
    molecular_orbitals: Vec<MolecularOrbital>,
}

pub struct MolecularOrbital {
    basis_functions: Vec<ContractedGaussian>,
    coefficients: DVector<f64>,
}

impl MolecularOrbital {
    pub fn new(basis_functions: Vec<ContractedGaussian>, coefficients: DVector<f64>) -> Self {
        Self {
            basis_functions,
            coefficients,
        }
    }

    pub fn evaluate(&self, at: &Vector3<f64>) -> f64 {
        self.coefficients
            .into_iter()
            .zip(self.basis_functions.iter())
            .map(|(coeff, basis)| coeff * basis.evaluate(at))
            .sum::<f64>()
    }
}

impl MolecularOrbitals {
    const MIN_COEFFICIENT_MAGNITUDE: f64 = 0.05;

    pub fn new(basis_functions: Vec<ContractedGaussian>, coeff_matrix: DMatrix<f64>) -> Self {
        let molecular_orbitals = coeff_matrix
            .column_iter()
            .map(|column| {
                let (indices, elements): (Vec<_>, Vec<_>) = column
                    .iter()
                    .enumerate()
                    .filter(|(_index, element)| element.abs() > Self::MIN_COEFFICIENT_MAGNITUDE)
                    .unzip();

                MolecularOrbital::new(
                    indices
                        .into_iter()
                        .map(|index| basis_functions[index].clone())
                        .collect(),
                    DVector::from_vec(elements),
                )
            })
            .collect();

        Self {
            coeff_matrix,
            molecular_orbitals,
        }
    }

    pub fn coeff_matrix(&self) -> &DMatrix<f64> {
        &self.coeff_matrix
    }
}

impl Index<usize> for MolecularOrbitals {
    type Output = MolecularOrbital;

    fn index(&self, index: usize) -> &Self::Output {
        &self.molecular_orbitals[index]
    }
}
