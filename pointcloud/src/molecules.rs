use basis_set::atom::Atom;
use basis_set::periodic_table::AtomType;
use basis_set::periodic_table::AtomType::{Hydrogen, Nitrogen, Oxygen};
use basis_set::BasisSet;
use nalgebra::Vector3;

#[allow(dead_code)]
/// angle: 104.5°
/// length: 1.8a_0
pub fn build_water(basis: &BasisSet, angle: f64, length: f64) -> [Atom; 3] {
    let a = angle.to_radians() * 0.5;
    [
        basis.get(Vector3::new(-a.sin(), -a.cos(), 0.0) * length, Hydrogen, 0),
        basis.get(Vector3::new(a.sin(), -a.cos(), 0.0) * length, Hydrogen, 0),
        basis.get(Vector3::zeros(), Oxygen, 0),
    ]
}

#[allow(dead_code)]
/// angle: 115°
/// length: 2.35
pub fn build_nitrite(basis: &BasisSet, angle: f64, length: f64) -> [Atom; 3] {
    let a = angle.to_radians() * 0.5;
    [
        basis.get(Vector3::new(-a.sin(), -a.cos(), 0.0) * length, Oxygen, 0),
        basis.get(Vector3::zeros() * length, Nitrogen, 1),
        basis.get(Vector3::new(a.sin(), -a.cos(), 0.0) * length, Oxygen, 0),
    ]
}

#[allow(dead_code)]
/// angle: 120°
/// length_cc: 2.5a_0
/// length_ch: 2.0a_0
pub fn build_ethene(basis: &BasisSet, angle: f64, length_cc: f64, length_ch: f64) -> [Atom; 6] {
    let a = (std::f64::consts::PI - angle.to_radians()) * 0.5;

    let carbon_a = basis.get(Vector3::new(length_cc * 0.5, 0.0, 0.0), AtomType::Carbon, 0);
    let carbon_b = basis.get(
        Vector3::new(-length_cc * 0.5, 0.0, 0.0),
        AtomType::Carbon,
        0,
    );

    [
        basis.get(
            carbon_a.position() + Vector3::new(a.sin(), 0.0, a.cos()) * length_ch,
            AtomType::Hydrogen,
            0,
        ),
        basis.get(
            carbon_a.position() + Vector3::new(a.sin(), 0.0, -a.cos()) * length_ch,
            AtomType::Hydrogen,
            0,
        ),
        carbon_a,
        basis.get(
            carbon_b.position() - Vector3::new(a.sin(), 0.0, a.cos()) * length_ch,
            AtomType::Hydrogen,
            0,
        ),
        basis.get(
            carbon_b.position() - Vector3::new(a.sin(), 0.0, -a.cos()) * length_ch,
            AtomType::Hydrogen,
            0,
        ),
        carbon_b,
    ]
}

#[allow(dead_code)]
/// length_cc: 2.25a_0
/// length_ch: 2.0a_0
pub fn build_ethyne(basis: &BasisSet, length_cc: f64, length_ch: f64) -> [Atom; 4] {
    [
        basis.get(
            Vector3::new(-length_ch * 0.5 - length_ch, 0.0, 0.0),
            AtomType::Hydrogen,
            0,
        ),
        basis.get(
            Vector3::new(-length_cc * 0.5, 0.0, 0.0),
            AtomType::Carbon,
            0,
        ),
        basis.get(Vector3::new(length_cc * 0.5, 0.0, 0.0), AtomType::Carbon, 0),
        basis.get(
            Vector3::new(length_ch * 0.5 + length_ch, 0.0, 0.0),
            AtomType::Hydrogen,
            0,
        ),
    ]
}

pub fn build_benzene(basis: &BasisSet, length_cc: f64, length_ch: f64) -> Vec<Atom> {
    let angle_incr = std::f64::consts::FRAC_PI_3;
    (0..6)
        .map(|k| k as f64 * angle_incr)
        .flat_map(|a| {
            [
                basis.get(
                    Vector3::new(a.sin(), 0.0, a.cos()) * length_cc,
                    AtomType::Carbon,
                    0,
                ),
                basis.get(
                    Vector3::new(a.sin(), 0.0, a.cos()) * (length_cc + length_ch),
                    AtomType::Hydrogen,
                    0,
                ),
            ]
        })
        .collect()
}
