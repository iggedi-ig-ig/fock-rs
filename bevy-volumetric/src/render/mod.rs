pub mod density;
pub mod molecule;
pub mod volume;

use bevy::prelude::*;

use self::{
    density::ComputeDensityPlugin, molecule::MoleculeRenderPlugin, volume::VolumeRenderPlugin,
};

#[derive(Resource)]
pub struct RenderSettings {
    pub current_energy_level: usize,
    pub density_resolution: usize,
    pub density_scale: f32,
}

impl Default for RenderSettings {
    fn default() -> Self {
        Self {
            current_energy_level: 0,
            density_resolution: 200,
            density_scale: 15.0,
        }
    }
}

pub struct RenderPlugin;

impl Plugin for RenderPlugin {
    fn build(&self, app: &mut App) {
        app.init_resource::<RenderSettings>().add_plugins((
            MoleculeRenderPlugin,
            ComputeDensityPlugin,
            VolumeRenderPlugin,
        ));
    }
}
