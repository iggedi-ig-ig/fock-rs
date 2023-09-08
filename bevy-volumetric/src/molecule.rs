use basis_set::atom::Atom;
use bevy::prelude::*;

#[derive(Resource, Deref, Default)]
pub struct LoadedMolecule(Option<Vec<Atom>>);

impl LoadedMolecule {
    pub fn load(&mut self, molecule: Vec<Atom>) {
        self.0 = Some(molecule);
    }
}

pub struct MoleculeLoaderPlugin;

impl Plugin for MoleculeLoaderPlugin {
    fn build(&self, app: &mut App) {
        app.insert_resource(LoadedMolecule::default());
    }
}
