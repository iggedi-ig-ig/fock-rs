pub mod contracted_gaussian;
pub mod primitives;
pub mod utils;

use nalgebra::Vector3;

pub trait BasisFunction {
    /// Overlap integral
    /// <a|b>
    fn overlap_int(a: &Self, b: &Self) -> f64;

    /// Kinetic energy integral
    /// <a|del^2|b>
    fn kinetic_int(a: &Self, b: &Self) -> f64;

    /// Nuclear attraction integral
    /// <a|1/eA|b>
    fn nuclear_attraction_int(
        a: &Self,
        b: &Self,
        nucleus_pos: &Vector3<f64>,
        nucleus_charge: f64,
    ) -> f64;

    /// Electron-Electron repulsion
    /// <ab|cd>
    fn electron_repulsion_int(a: &Self, b: &Self, c: &Self, d: &Self) -> f64;

    fn evaluate(&self, at: &Vector3<f64>) -> f64;
}
