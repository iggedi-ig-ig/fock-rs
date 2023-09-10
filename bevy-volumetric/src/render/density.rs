use std::{ops::Deref, time::Instant};

use bevy::{
    prelude::*,
    render::{
        extract_resource::{ExtractResource, ExtractResourcePlugin},
        render_resource::{Extent3d, TextureDimension, TextureFormat, TextureUsages},
    },
    tasks::{AsyncComputeTaskPool, Task},
};
use futures_lite::future;
use nalgebra::Vector3;

use crate::hf::ConvergedScf;

use super::RenderSettings;

#[derive(Resource, Default)]
pub struct DensityBuffer(Vec<Option<Vec<f32>>>);

impl DensityBuffer {
    pub fn allocate(size: usize) -> Self {
        Self(vec![None; size])
    }

    // pub fn allocate_level(&mut self, n: usize) -> &mut Vec<f32> {
    //     let capacity = n.pow(3);
    //     self.0[n] = Some(vec![0.0; capacity]);
    //     self.0[n].as_mut().unwrap()
    // }

    pub fn level_mut(&mut self, n: usize) -> &mut Option<Vec<f32>> {
        &mut self.0[n]
    }

    pub fn level(&self, n: usize) -> Option<&[f32]> {
        self.0[n].as_deref()
    }
}

pub struct ComputeDensityPlugin;

impl Plugin for ComputeDensityPlugin {
    fn build(&self, app: &mut App) {
        app.init_resource::<DensityBuffer>()
            .init_resource::<Density3dTexture>()
            .add_plugins(ExtractResourcePlugin::<Density3dTexture>::default())
            .add_systems(
                Update,
                (update_densities, handle_density_tasks, update_3d_texture),
            );
    }
}

#[derive(Component)]
struct UpdateDensityTask {
    energy_level: usize,
    start: Instant,
    task: Task<Vec<f32>>,
}

fn update_densities(
    mut commands: Commands,
    mut density_buffer: ResMut<DensityBuffer>,
    converged: Res<ConvergedScf>,
    settings: Query<&RenderSettings>,
) {
    let settings = settings.single();
    let pool = AsyncComputeTaskPool::get();

    if converged.is_changed() {
        let ConvergedScf::Solved(result) = converged.deref() else {
            return;
        };

        *density_buffer = DensityBuffer::allocate(result.n_basis);
        let res = settings.density_resolution as usize;

        for energy_level in 0..result.n_basis {
            let orbitals = result.orbitals.clone();
            commands.spawn(UpdateDensityTask {
                energy_level,
                start: Instant::now(),
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
                "computed density buffer for energy level {}. Took {:?}",
                update_density_task.energy_level,
                update_density_task.start.elapsed()
            );

            commands.entity(entity).despawn_recursive();
        }
    }
}

#[derive(Resource, ExtractResource, Deref, Default, DerefMut, Clone)]
pub struct Density3dTexture(Handle<Image>);

fn update_3d_texture(
    mut density_texture: ResMut<Density3dTexture>,
    mut images: ResMut<Assets<Image>>,
    render_settings: Query<&RenderSettings>,
    density_buffer: Res<DensityBuffer>,
) {
    if density_buffer.is_changed() {
        let render_settings = render_settings.single();
        if images.contains(&density_texture.0) {
            images.remove(density_texture.0.clone());
        }

        let mut image = Image::new_fill(
            Extent3d {
                width: render_settings.density_resolution,
                height: render_settings.density_resolution,
                depth_or_array_layers: render_settings.density_resolution,
            },
            TextureDimension::D3,
            &(0.0_f32).to_be_bytes(),
            TextureFormat::R32Float,
        );
        image.texture_descriptor.usage = TextureUsages::COPY_DST | TextureUsages::STORAGE_BINDING;
        **density_texture = images.add(image);
    }
}
