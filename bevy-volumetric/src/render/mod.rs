pub mod density;
pub mod molecule;
pub mod volume;

use bevy::prelude::*;

use crate::hf::ConvergedScf;

use self::{
    density::ComputeDensityPlugin, molecule::MoleculeRenderPlugin, volume::VolumeRenderPlugin,
};

#[derive(Resource, Clone, Copy, Default)]
pub struct RenderSettings {
    pub current_energy_level: u32,
}

pub struct RenderPlugin;

impl Plugin for RenderPlugin {
    fn build(&self, app: &mut App) {
        app.add_plugins((
            MoleculeRenderPlugin,
            ComputeDensityPlugin,
            VolumeRenderPlugin,
        ))
        .add_systems(Update, update_energy_level)
        .init_resource::<RenderSettings>();
    }
}

fn update_energy_level(
    mut render_settings: ResMut<RenderSettings>,
    keys: Res<Input<KeyCode>>,
    converged: Res<ConvergedScf>,
) {
    if let ConvergedScf::Solved(sol) = &*converged {
        if keys.just_pressed(KeyCode::Left) {
            render_settings.current_energy_level = render_settings
                .current_energy_level
                .checked_sub(1)
                .unwrap_or(sol.n_basis as _);
        }
        if keys.just_pressed(KeyCode::Right) {
            render_settings.current_energy_level =
                (render_settings.current_energy_level + 1) % sol.n_basis as u32;
        }
    }
}
