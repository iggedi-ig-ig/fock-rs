pub mod contracted_gaussian;
pub mod electron_tensor;
pub mod primitives;
pub mod utils;

use nalgebra::Vector3;

#[derive(Copy, Clone, Debug)]
pub struct PointCharge {
    pub position: Vector3<f64>,
    pub charge: f64,
}

/// Trait for a basis function used in quantum chemistry calculations.
pub trait BasisFunction {
    /// Calculates the overlap integral <a|b>, where a and b are two basis functions.
    ///
    /// The overlap integral represents the amount of overlap between the two basis functions a and b.
    /// The value of the integral ranges from 0 to 1, with 1 indicating complete overlap.
    fn overlap_int(a: &Self, b: &Self) -> f64;

    /// Calculates the kinetic energy integral <a|del^2|b>, where a and b are two basis functions.
    ///
    /// The kinetic energy integral represents the kinetic energy of an electron moving between the two
    /// basis functions a and b. The value of the integral is always positive.
    fn kinetic_int(a: &Self, b: &Self) -> f64;

    /// Calculates the nuclear attraction integral <a|1/eA|b>, where a and b are two basis functions,
    /// and eA is the electric charge of the atomic nucleus A.
    ///
    /// The nuclear attraction integral represents the electrostatic attraction between the electron
    /// in basis function b and the atomic nucleus A. The value of the integral is always negative.
    /// This is because electrons are attracted to the positively charged atomic nucleus.
    fn nuclear_attraction_int(a: &Self, b: &Self, nuclei: &[PointCharge]) -> f64;

    /// Calculates the electron-electron repulsion integral <ab|cd>, where a, b, c, and d are four basis functions.
    ///
    /// The electron-electron repulsion integral represents the electrostatic repulsion between two electrons
    /// in basis functions a and b, and two electrons in basis functions c and d. The value of the integral
    /// is always positive.
    fn electron_repulsion_int(a: &Self, b: &Self, c: &Self, d: &Self) -> f64;

    /// Evaluates the basis function at a given point in space.
    ///
    /// Given a basis function, this method returns the value of the basis function at a given point in space.
    fn evaluate(&self, at: &Vector3<f64>) -> f64;
}
