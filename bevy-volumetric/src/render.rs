use futures_lite::future;
use nalgebra::Vector3;
use std::ops::Deref;

use bevy::{
    prelude::*,
    tasks::{AsyncComputeTaskPool, Task},
};

use crate::{hf::ConvergedScf, molecule::LoadedMolecule};

#[derive(Resource, Default)]
pub struct DensityBuffer(Vec<Option<Vec<f32>>>);

impl DensityBuffer {
    pub fn allocate(size: usize) -> Self {
        Self(vec![None; size])
    }

    pub fn allocate_level(&mut self, n: usize) -> &mut Vec<f32> {
        let capacity = n.pow(3);
        self.0[n] = Some(vec![0.0; capacity]);
        self.0[n].as_mut().unwrap()
    }

    pub fn level_mut(&mut self, n: usize) -> &mut Option<Vec<f32>> {
        &mut self.0[n]
    }
}

#[derive(Resource)]
pub struct RenderSettings {
    pub density_resolution: usize,
    pub density_scale: f32,
}

impl Default for RenderSettings {
    fn default() -> Self {
        Self {
            density_resolution: 200,
            density_scale: 15.0,
        }
    }
}

pub struct RenderPlugin;

impl Plugin for RenderPlugin {
    fn build(&self, app: &mut App) {
        app.init_resource::<RenderSettings>()
            .init_resource::<DensityBuffer>()
            .add_systems(
                Update,
                (update_densities, handle_density_tasks, add_molecule_atoms),
            );
    }
}

#[derive(Component)]
struct AtomMarker;

fn add_molecule_atoms(
    mut commands: Commands,
    mut meshes: ResMut<Assets<Mesh>>,
    mut materials: ResMut<Assets<StandardMaterial>>,
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
        let [red, green, blue] = atom.atom_type().color();
        commands.spawn((
            PbrBundle {
                mesh: meshes.add(Mesh::from(shape::UVSphere::default())),
                material: materials.add(StandardMaterial {
                    base_color: Color::Rgba {
                        red,
                        green,
                        blue,
                        alpha: 1.0,
                    },
                    ..Default::default()
                }),
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

#[derive(Component)]
struct UpdateDensityTask {
    energy_level: usize,
    task: Task<Vec<f32>>,
}

fn update_densities(
    mut commands: Commands,
    mut density_buffer: ResMut<DensityBuffer>,
    converged: Res<ConvergedScf>,
    settings: Res<RenderSettings>,
) {
    let pool = AsyncComputeTaskPool::get();

    if converged.is_changed() {
        let ConvergedScf::Solved(result) = converged.deref() else {
           return;
        };

        *density_buffer = DensityBuffer::allocate(result.n_basis);
        let res = settings.density_resolution;

        for energy_level in 0..result.n_basis {
            let orbitals = result.orbitals.clone();
            commands.spawn(UpdateDensityTask {
                energy_level,
                task: pool.spawn(async move {
                    let mut output = vec![0.0; res.pow(3)];
                    for (i, entry) in output.iter_mut().enumerate() {
                        let x = i % res;
                        let y = (i / res) % res;
                        let z = i / res.pow(2);

                        let at = Vector3::new(x, y, z).map(|x| x as f64 / res as f64 - 0.5);
                        *entry = orbitals[energy_level].evaluate(&at) as f32;
                    }
                    output
                }),
            });
        }
    }
}

fn handle_density_tasks(
    mut commands: Commands,
    mut query: Query<(Entity, &mut UpdateDensityTask)>,
    mut density_buffer: ResMut<DensityBuffer>,
) {
    for (entity, mut update_density_task) in &mut query {
        if let Some(buffer) = future::block_on(future::poll_once(&mut update_density_task.task)) {
            *density_buffer.level_mut(update_density_task.energy_level) = Some(buffer);
            info!(
                "computed density buffer for energy level {}",
                update_density_task.energy_level
            );

            commands.entity(entity).despawn_recursive();
        }
    }
}
