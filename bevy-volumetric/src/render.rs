use std::ops::Deref;
use nalgebra::Vector3;

use bevy::prelude::*;

use crate::scf::ConvergedScf;

#[derive(Resource)]
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
        app.insert_resource(RenderSettings::default())
            .add_systems(Update, update_densities);
    }
}

fn update_densities(
    mut density_texture: ResMut<DensityBuffer>,
    converged: Res<ConvergedScf>,
    settings: Res<RenderSettings>
) {
    // TODO: choose energy level
    let energy_level= 0;
    if converged.is_changed() {
        let ConvergedScf::Solved(result) = converged.deref() else {
           return; 
        };

        let res = settings.density_resolution;
        let output = density_texture.allocate_level(energy_level);
        for (i, entry) in output.iter_mut().enumerate() {
            let x = i % res;
            let y = (i / res) % res;
            let z = i / res.pow(2);

            let at = Vector3::new(x, y, z).map(|x| x as f64 / res as f64 - 0.5);
            *entry = result.orbitals[energy_level].evaluate(&at) as f32;
        }
    }
}
