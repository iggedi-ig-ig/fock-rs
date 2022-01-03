#![allow(dead_code)]

pub mod basis_sets;

use basis::contracted_gaussian::ContractedGaussian;
use basis::primitives::GaussianPrimitive;
use nalgebra::Vector3;
use serde::Deserialize;
use std::collections::HashMap;
use std::fmt::{Display, Formatter};

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
    elements: HashMap<String, ElectronicConfiguration>,
    version: String,
    function_types: Vec<FunctionType>,
    names: Vec<String>,
    tags: Vec<String>,
    family: BasisSetFamily,
    description: String,
    role: String,
    name: String,
}

#[derive(Clone, Debug)]
pub struct Atom {
    position: Vector3<f64>,
    ordinal: usize,
    oxidation_state: i32,
    basis: Vec<ContractedGaussian>,
}

impl Atom {
    pub fn new(
        position: Vector3<f64>,
        ordinal: usize,
        basis: Vec<ContractedGaussian>,
        oxidation_state: i32,
    ) -> Self {
        Self {
            position,
            ordinal,
            oxidation_state,
            basis,
        }
    }

    pub fn position(&self) -> Vector3<f64> {
        self.position
    }
    pub fn ordinal(&self) -> usize {
        self.ordinal
    }
    pub fn basis(&self) -> &[ContractedGaussian] {
        &self.basis
    }
    pub fn num_electrons(&self) -> usize {
        (self.ordinal as i32 + self.oxidation_state) as usize
    }
}

impl BasisSet {
    pub fn get(
        &self,
        position: Vector3<f64>,
        ordinal: usize,
        oxidation_state: i32,
    ) -> Option<Atom> {
        let config = self.elements.get(&ordinal.to_string());
        config.map(|config| {
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
            Atom::new(position, ordinal, basis_functions, oxidation_state)
        })
    }
}

impl Display for Atom {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        writeln!(
            f,
            "Atom Num {} (electrons: {})\n\tPos: {:0.3}, {:0.3}, {:0.3}\n\tBasis Functions:",
            self.ordinal,
            self.num_electrons(),
            self.position.x,
            self.position.y,
            self.position.z
        )?;
        for (i, contracted_gaussian) in self.basis.iter().enumerate() {
            writeln!(f, "\t{}:", i)?;
            for primitive in contracted_gaussian.primitives() {
                writeln!(f, "\t\t{}", primitive)?;
            }
        }
        Ok(())
    }
}

#[test]
pub fn test_basis() {
    use crate::basis_sets::BASIS_6_31G;

    if let Some(atom) = BASIS_6_31G.get(Vector3::zeros(), 6, 0) {
        println!("{}", atom)
    } else {
        panic!("Couldn't get hydrogen");
    }
}
