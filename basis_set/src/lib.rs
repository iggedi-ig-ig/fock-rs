#![allow(dead_code)]

pub mod atom;
pub mod basis_sets;
pub mod periodic_table;

use crate::atom::Atom;
use crate::periodic_table::AtomType;
use basis::contracted_gaussian::ContractedGaussian;
use basis::primitives::GaussianPrimitive;
use nalgebra::Vector3;
use serde::Deserialize;
use std::collections::HashMap;

#[derive(Deserialize, Debug)]
pub struct MolssiBseSchema {
    schema_type: String,
    schema_version: String,
}

#[derive(Deserialize, Debug)]
pub struct ElectronShell {
    function_type: String,
    region: String,
    angular_momentum: Vec<i32>,
    exponents: Vec<String>,
    coefficients: Vec<Vec<String>>,
}

#[derive(Deserialize, Debug)]
struct Reference {
    reference_description: String,
    reference_keys: Vec<String>,
}

#[derive(Deserialize, Debug)]
pub struct ElectronicConfiguration {
    electron_shells: Vec<ElectronShell>,
    references: Vec<Reference>,
}

#[derive(Deserialize, Debug)]
pub enum FunctionType {
    #[serde(rename = "gto")]
    Gaussian,
    #[serde(rename = "gto_spherical")]
    GaussianSpherical,
    #[serde(rename = "gto_cartesian")]
    GaussianCartesian,
}

#[derive(Deserialize, Debug)]
pub enum BasisSetFamily {
    #[serde(rename = "pople")]
    Pople,
    #[serde(rename = "sto")]
    SlaterType,
}

#[derive(Deserialize, Debug)]
pub struct BasisSet {
    molssi_bse_schema: MolssiBseSchema,
    revision_description: String,
    revision_date: String,
    elements: HashMap<AtomType, ElectronicConfiguration>,
    version: String,
    function_types: Vec<FunctionType>,
    names: Vec<String>,
    tags: Vec<String>,
    family: BasisSetFamily,
    description: String,
    role: String,
    name: String,
}

impl BasisSet {
    pub fn get(&self, position: Vector3<f64>, atom_type: AtomType) -> Atom {
        self.elements
            .get(&atom_type)
            .map(|config| {
                let mut basis_functions = Vec::new();
                for shell in config.electron_shells.iter() {
                    for (idx, angular) in shell.angular_momentum.iter().enumerate() {
                        let angulars = (0..=*angular)
                            .flat_map(|x| {
                                (0..=*angular).flat_map(move |y| {
                                    (0..=*angular).filter_map(move |z| {
                                        if x + y + z == *angular {
                                            Some([x, y, z])
                                        } else {
                                            None
                                        }
                                    })
                                })
                            })
                            .collect::<Vec<[i32; 3]>>();

                        for angular in angulars {
                            basis_functions.push(ContractedGaussian::new(
                                position,
                                shell
                                    .exponents
                                    .iter()
                                    .zip(shell.coefficients[idx].iter())
                                    .map(|(s1, s2)| [s1.parse(), s2.parse()].map(Result::unwrap))
                                    .map(|[exp, coeff]| GaussianPrimitive::new(angular, exp, coeff))
                                    .collect(),
                            ))
                        }
                    }
                }
                Atom::new(position, atom_type, basis_functions)
            })
            .expect("Failed to get atom from basis set")
    }
}

#[test]
pub fn test_basis() {
    use crate::basis_sets::BASIS_6_31G;

    println!(
        "{}",
        BASIS_6_31G.get(Vector3::zeros(), AtomType::Phosphorous)
    )
}
