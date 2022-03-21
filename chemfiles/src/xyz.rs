use basis_set::atom::Atom;
use basis_set::periodic_table::AtomType;
use basis_set::BasisSet;
use nalgebra::Vector3;
use regex::Regex;
use std::fs::File;
use std::io::Read;
use std::path::Path;

/// XYZ File format

pub fn read_xyz_file<P: AsRef<Path>>(
    path: P,
    basis_set: &BasisSet,
) -> Result<Vec<Atom>, Box<dyn std::error::Error>> {
    let pattern = Regex::new(r"\s*(\w+)\s+([+\-0-9.]+)\s+([+\-0-9.]+)\s+([0-9.+\-]+).*\b")?;

    let mut file = File::open(path)?;
    let mut content = String::new();

    file.read_to_string(&mut content)?;

    let mut atoms = Vec::new();
    for cap in pattern.captures_iter(&content) {
        println!("{cap:?}");

        let symbol = cap[1].to_uppercase();
        let x = cap[2].parse::<f64>()?;
        let y = cap[3].parse::<f64>()?;
        let z = cap[4].parse::<f64>()?;

        println!("{symbol} {x} {y} {z}");
        atoms.push(basis_set.get(Vector3::new(x, y, z), AtomType::from_symbol(&symbol)));
    }

    Ok(atoms)
}