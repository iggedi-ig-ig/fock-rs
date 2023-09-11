use bevy::{
    pbr::{NotShadowCaster, NotShadowReceiver},
    prelude::*,
    reflect::{TypePath, TypeUuid},
    render::render_resource::{AsBindGroup, ShaderType},
};

use super::density::Density3dTexture;

pub struct VolumeRenderPlugin;

impl Plugin for VolumeRenderPlugin {
    fn build(&self, app: &mut App) {
        app.add_plugins(MaterialPlugin::<VolumetricMaterial> {
            prepass_enabled: false,
            ..Default::default()
        })
        .add_systems(Startup, setup)
        .add_systems(Update, update_density_texture)
        .init_resource::<VolumetricSettings>();
    }
}

fn setup(
    mut commands: Commands,
    mut meshes: ResMut<Assets<Mesh>>,
    mut materials: ResMut<Assets<VolumetricMaterial>>,
    settings: Res<VolumetricSettings>,
) {
    let scale = settings.box_max - settings.box_min;

    commands.spawn((
        MaterialMeshBundle {
            mesh: meshes.add(Mesh::from(shape::Cube::default())),
            material: materials.add(VolumetricMaterial {
                settings: *settings,
                density_texture: None,
            }),
            transform: Transform::from_scale(scale),
            ..Default::default()
        },
        NotShadowCaster,
        NotShadowReceiver,
    ));
}

fn update_density_texture(
    mut materials: ResMut<Assets<VolumetricMaterial>>,
    density_texture: Res<Density3dTexture>,
    volume_materials: Query<&Handle<VolumetricMaterial>>,
) {
    if !density_texture.is_changed() {
        return;
    }

    for volume_material in &volume_materials {
        let Some(material) = materials.get_mut(volume_material) else {
            continue;
        };

        info!("update density texture");
        material.density_texture = Some((**density_texture).clone());
    }
}

#[derive(Resource, Debug, Copy, Clone, ShaderType)]
pub struct VolumetricSettings {
    pub resolution: u32,
    pub box_min: Vec3,
    pub box_max: Vec3,
}

impl VolumetricSettings {
    pub const SCALE: f32 = 10.0;
}

impl Default for VolumetricSettings {
    fn default() -> Self {
        Self {
            resolution: 100,
            box_min: -Vec3::ONE * Self::SCALE * 0.5,
            box_max: Vec3::ONE * Self::SCALE * 0.5,
        }
    }
}

#[derive(TypePath, TypeUuid, AsBindGroup, Debug, Clone)]
#[uuid = "3b3dfb6b-2683-4b7f-978d-77d7d9acde34"]
pub struct VolumetricMaterial {
    #[uniform(0)]
    settings: VolumetricSettings,
    #[texture(1, dimension = "3d")]
    #[sampler(2)]
    density_texture: Option<Handle<Image>>,
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
