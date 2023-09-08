use basis_set::atom::Atom;
use bevy::prelude::*;

#[derive(Resource, Deref, Default)]
pub struct LoadedMolecule(Option<Vec<Atom>>);

pub struct MoleculeLoaderPlugin;

impl Plugin for MoleculeLoaderPlugin {
    fn build(&self, app: &mut App) {
        app.insert_resource(LoadedMolecule::default());
    }
}
