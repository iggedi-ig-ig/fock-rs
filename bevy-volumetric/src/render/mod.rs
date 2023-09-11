pub mod density;
pub mod molecule;
pub mod volume;

use bevy::prelude::*;

use self::{
    density::ComputeDensityPlugin, molecule::MoleculeRenderPlugin, volume::VolumeRenderPlugin,
};

#[derive(Resource, Clone, Copy)]
pub struct RenderSettings {
    pub current_energy_level: u32,
}

impl Default for RenderSettings {
    fn default() -> Self {
        Self {
            current_energy_level: 2,
        }
    }
}

pub struct RenderPlugin;

impl Plugin for RenderPlugin {
    fn build(&self, app: &mut App) {
        app.add_plugins((
            MoleculeRenderPlugin,
            ComputeDensityPlugin,
            VolumeRenderPlugin,
        ))
        .init_resource::<RenderSettings>();
    }
}
