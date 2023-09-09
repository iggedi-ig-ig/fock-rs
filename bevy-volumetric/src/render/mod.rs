mod density;
mod molecule;

use bevy::prelude::*;

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
        app.init_resource::<RenderSettings>();
    }
}
