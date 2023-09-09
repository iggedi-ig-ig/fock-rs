use basis_set::periodic_table::AtomType;
use bevy::{prelude::*, utils::HashMap};

use crate::molecule::LoadedMolecule;

pub struct MoleculeRenderPlugin;

impl Plugin for MoleculeRenderPlugin {
    fn build(&self, app: &mut App) {
        app.add_systems(Startup, setup)
            .add_systems(Update, add_molecule_atoms);
    }
}

#[derive(Resource, Deref)]
struct AtomMesh(Handle<Mesh>);

#[derive(Resource, Deref, DerefMut, Default)]
struct AtomMaterial(HashMap<AtomType, Handle<StandardMaterial>>);

fn setup(mut commands: Commands, mut meshes: ResMut<Assets<Mesh>>) {
    let atom_mesh = Mesh::from(shape::UVSphere::default());
    commands.insert_resource(AtomMesh(meshes.add(atom_mesh)));
    commands.init_resource::<AtomMaterial>();
}

#[derive(Component)]
struct AtomMarker;

fn add_molecule_atoms(
    mut commands: Commands,
    mut materials: ResMut<Assets<StandardMaterial>>,
    mut atom_material: ResMut<AtomMaterial>,
    atom_mesh: Res<AtomMesh>,
    old_query: Query<Entity, With<AtomMarker>>,
    loaded_molecule: Res<LoadedMolecule>,
) {
    if !loaded_molecule.is_changed() {
        return;
    }

    for entity in &old_query {
        commands.entity(entity).despawn_recursive();
    }

    let Some(loaded_molecule) = loaded_molecule.as_deref() else {
        return;
    };

    for atom in loaded_molecule {
        let mesh = (*atom_mesh).clone();
        let material = atom_material
            .entry(atom.atom_type())
            .or_insert_with(|| {
                let [red, green, blue] = atom.atom_type().color();

                materials.add(StandardMaterial {
                    base_color: Color::Rgba {
                        red,
                        green,
                        blue,
                        alpha: 1.0,
                    },
                    ..Default::default()
                })
            })
            .clone();
        commands.spawn((
            PbrBundle {
                mesh,
                material,
                transform: Transform::from_xyz(
                    atom.position().x as _,
                    atom.position().y as _,
                    atom.position().z as _,
                )
                .with_scale(Vec3::ONE * 0.5),
                ..Default::default()
            },
            AtomMarker,
        ));
    }
}
