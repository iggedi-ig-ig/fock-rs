use bevy::{
    pbr::{NotShadowCaster, NotShadowReceiver},
    prelude::*,
    reflect::{TypePath, TypeUuid},
    render::render_resource::{AsBindGroup, ShaderType},
};

pub struct VolumeRenderPlugin;

impl Plugin for VolumeRenderPlugin {
    fn build(&self, app: &mut App) {
        app.add_plugins(MaterialPlugin::<VolumetricMaterial> {
            prepass_enabled: false,
            ..Default::default()
        })
        .add_systems(Startup, setup);
    }
}

fn setup(
    mut commands: Commands,
    mut meshes: ResMut<Assets<Mesh>>,
    mut materials: ResMut<Assets<VolumetricMaterial>>,
) {
    commands.spawn((
        MaterialMeshBundle {
            mesh: meshes.add(Mesh::from(shape::Cube::default())),
            material: materials.add(VolumetricMaterial {
                settings: VolumetricSettings {
                    resolution: 100,
                    box_size: 10.0,
                },
            }),
            transform: Transform::IDENTITY
                .with_translation(-Vec3::ONE * 0.5)
                .with_scale(Vec3::ONE * 10.0),
            ..Default::default()
        },
        NotShadowCaster,
        NotShadowReceiver,
    ));
}

#[derive(Debug, Clone, Default, ShaderType)]
pub struct VolumetricSettings {
    resolution: u32,
    box_size: f32,
}

#[derive(TypePath, TypeUuid, AsBindGroup, Debug, Clone)]
#[uuid = "3b3dfb6b-2683-4b7f-978d-77d7d9acde34"]
pub struct VolumetricMaterial {
    #[uniform(0)]
    settings: VolumetricSettings,
}

impl Material for VolumetricMaterial {
    fn fragment_shader() -> bevy::render::render_resource::ShaderRef {
        "shaders/volumetric.wgsl".into()
    }

    fn alpha_mode(&self) -> AlphaMode {
        AlphaMode::Blend
    }

    fn specialize(
        _pipeline: &bevy::pbr::MaterialPipeline<Self>,
        descriptor: &mut bevy::render::render_resource::RenderPipelineDescriptor,
        _layout: &bevy::render::mesh::MeshVertexBufferLayout,
        _key: bevy::pbr::MaterialPipelineKey<Self>,
    ) -> Result<(), bevy::render::render_resource::SpecializedMeshPipelineError> {
        descriptor.primitive.cull_mode = None;
        Ok(())
    }
}
