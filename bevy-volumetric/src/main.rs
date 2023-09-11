pub mod hf;
pub mod molecule;
pub mod render;

use bevy::{core_pipeline::prepass::DepthPrepass, prelude::*};
use hf::HfPlugin;
use molecule::{LoadedMolecule, MoleculeLoaderPlugin};
use render::RenderPlugin;

fn main() {
    App::new()
        .add_plugins((DefaultPlugins, MoleculeLoaderPlugin, RenderPlugin, HfPlugin))
        .add_systems(Startup, setup)
        .add_systems(Update, update_cam_pos)
        .insert_resource(ClearColor(Color::GRAY))
        .run();
}

fn setup(mut commands: Commands, mut molecule: ResMut<LoadedMolecule>) {
    commands.spawn((
        Camera3dBundle {
            transform: Transform::from_xyz(0.0, 10.0, -20.0).looking_at(Vec3::ZERO, Vec3::Y),
            ..Default::default()
        },
        DepthPrepass,
    ));
    commands.spawn(PointLightBundle {
        point_light: PointLight {
            color: Color::WHITE,
            intensity: 1000.0,
            ..Default::default()
        },
        transform: Transform::from_xyz(0.0, 5.0, 0.0),
        ..Default::default()
    });

    molecule.load(
        chemfiles::xyz::read_xyz_file(
            "chemfiles/molecules/water.xyz",
            &basis_set::basis_sets::BASIS_STO_3G,
        )
        .expect("failed to read molecule"),
    );
}

fn update_cam_pos(mut camera: Query<&mut Transform, With<Camera>>, time: Res<Time>) {
    let mut camera = camera.single_mut();

    camera.rotate_around(
        Vec3::ZERO,
        Quat::from_axis_angle(Vec3::Y, 0.25 * time.delta_seconds()),
    );
}
