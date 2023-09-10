pub mod density;
pub mod molecule;
pub mod volume;

use bevy::{
    prelude::*,
    render::{extract_component::ExtractComponent, render_resource::ShaderType},
};

use self::{
    density::ComputeDensityPlugin, molecule::MoleculeRenderPlugin, volume::VolumeRenderPlugin,
};

#[derive(Component, Clone, Copy, ExtractComponent, ShaderType)]
pub struct RenderSettings {
    pub current_energy_level: u32,
    pub density_resolution: u32,
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
        app.add_plugins((
            MoleculeRenderPlugin,
            ComputeDensityPlugin,
            VolumeRenderPlugin,
        ))
        .add_systems(Startup, setup);
    }
}

fn setup(mut commands: Commands) {
    commands.spawn(RenderSettings::default());
}
