pub mod hf;
pub mod molecule;
pub mod render;

use bevy::prelude::*;
use hf::HfPlugin;
use molecule::{LoadedMolecule, MoleculeLoaderPlugin};
use render::RenderPlugin;

fn main() {
    App::new()
        .add_plugins((DefaultPlugins, MoleculeLoaderPlugin, RenderPlugin, HfPlugin))
        .add_systems(Startup, setup)
        .insert_resource(ClearColor(Color::GRAY))
        .run();
}

fn setup(mut commands: Commands, mut molecule: ResMut<LoadedMolecule>) {
    commands.spawn(Camera3dBundle::default());

    molecule.load(
        chemfiles::xyz::read_xyz_file(
            "chemfiles/molecules/benzene.xyz",
            &basis_set::basis_sets::BASIS_6_31G,
        )
        .expect("failed to read molecule"),
    );
}
