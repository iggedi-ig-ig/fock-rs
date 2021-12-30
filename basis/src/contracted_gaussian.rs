use crate::primitives::GaussianPrimitive;
use crate::BasisFunction;
use nalgebra::Vector3;

pub struct ContractedGaussian<const N: usize> {
    position: Vector3<f64>,
    primitives: [GaussianPrimitive; N],
}

impl<const N: usize> ContractedGaussian<N> {
    pub fn position(&self) -> Vector3<f64> {
        self.position
    }
    pub fn primitives(&self) -> &[GaussianPrimitive; N] {
        &self.primitives
    }
}

impl<const N: usize> BasisFunction for ContractedGaussian<N> {
    fn overlap_int(a: &Self, b: &Self) -> f64 {
        todo!("Overlap integral isn't implemented for ContractedGaussians yet")
    }

    fn kinetic_int(a: &Self, b: &Self) -> f64 {
        todo!("Kinetic energy integral isn't implemented for ContractedGaussians yet")
    }

    fn nuclear_attraction_int(
        a: &Self,
        b: &Self,
        nucleus_pos: &Vector3<f64>,
        nucleus_charge: f64,
    ) -> f64 {
        todo!("Nuclear attraction energy integral isn't implemented for ContractedGaussians yet")
    }

    fn electron_repulsion_int(a: &Self, b: &Self, c: &Self, d: &Self) -> f64 {
        todo!("Electron-Electron repulsion energy integral isn't implemented for ContractedGaussians yet")
    }
}
