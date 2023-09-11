#import bevy_pbr::mesh_types
#import bevy_pbr::mesh_view_bindings  globals
#import bevy_pbr::prepass_utils
#import bevy_pbr::mesh_vertex_output  MeshVertexOutput
#import bevy_pbr::pbr_functions as pbr_functions

struct VolumetricSettings {
    resolution: u32,
    box_size: f32
}
@group(1) @binding(0)
var<uniform> settings: VolumetricSettings;

@fragment
fn fragment(
    mesh: MeshVertexOutput,
) -> @location(0) vec4<f32> {
    let depth = bevy_pbr::prepass_utils::prepass_depth(mesh.position, 0u);
    let v = pbr_functions::calculate_view(mesh.world_position, false);
    let start = v * depth;
    return vec4(vec3(depth), 1.0);
}