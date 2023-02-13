use basis::contracted_gaussian::ContractedGaussian;
use basis::BasisFunction;
use nalgebra::{DMatrix, DVector, Vector3};
use std::ops::Index;

pub struct MolecularOrbitals {
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

    /// Evaluates the function defined by this struct at the given point.
    ///
    /// Computes the value of the function by multiplying each coefficient with the
    /// corresponding basis function evaluated at the given point `at`. The resulting
    /// sequence of values is summed to produce the final output value of the function.
    ///
    /// # Arguments
    ///
    /// * `at` - A reference to a `Vector3<f64>` object representing the point at which
    ///          to evaluate the function.
    ///
    /// # Returns
    ///
    /// The value of the function at the given point, as a `f64`.
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

    /// Creates a new set of molecular orbitals from the given basis functions and coefficient matrix.
    ///
    /// The method constructs a set of `MolecularOrbital` objects from the given `basis_functions`
    /// and the columns of the `coeff_matrix` that have coefficients with magnitudes greater than
    /// `MIN_COEFFICIENT_MAGNITUDE`. Each `MolecularOrbital` object is created by selecting the
    /// relevant basis functions and coefficients from the `basis_functions` and `coeff_matrix`
    /// inputs, respectively. The resulting set of `MolecularOrbital` objects is returned as a new
    /// `MolecularOrbitals` struct.
    ///
    /// # Arguments
    ///
    /// * `basis_functions` - A vector of `ContractedGaussian` objects representing the basis functions
    ///                       used to construct the molecular orbitals.
    /// * `coeff_matrix` - A reference to a `DMatrix<f64>` object representing the coefficient matrix
    ///                    used to weight the basis functions when constructing the molecular orbitals.
    ///
    /// # Returns
    ///
    /// A new `MolecularOrbitals` struct containing the set of molecular orbitals constructed from the
    /// input basis functions and coefficient matrix.
    pub fn new(basis_functions: Vec<ContractedGaussian>, coeff_matrix: &DMatrix<f64>) -> Self {
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

        Self { molecular_orbitals }
    }
}

impl Index<usize> for MolecularOrbitals {
    type Output = MolecularOrbital;

    fn index(&self, index: usize) -> &Self::Output {
        &self.molecular_orbitals[index]
    }
}
