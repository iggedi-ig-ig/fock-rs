#import bevy_pbr::mesh_types
#import bevy_pbr::mesh_view_bindings view
#import bevy_pbr::prepass_utils
#import bevy_pbr::mesh_vertex_output  MeshVertexOutput
#import bevy_pbr::pbr_functions as pbr_functions

struct VolumetricSettings {
    resolution: u32,
    density_multiplier: f32,
    box_size: f32,
    box_min: vec3<f32>,
    box_max: vec3<f32>,
}
@group(1) @binding(0)
var<uniform> settings: VolumetricSettings;
@group(1) @binding(1)
var density_texture: texture_3d<f32>;
@group(1) @binding(2)
var density_sampler: sampler;

// ported from hlsl of Sebastian Lague's Cloud project
fn ray_box(ray_origin: vec3<f32>, dir: vec3<f32>, bounds_min: vec3<f32>, bounds_max: vec3<f32>) -> vec2<f32> {
    let inv_dir = 1.0 / dir;
    let t0 = (bounds_min - ray_origin) * inv_dir;
    let t1 = (bounds_max - ray_origin) * inv_dir;

    let tmin = min(t0, t1);
    let tmax = max(t0, t1);

    let dst_a = max(max(tmin.x, tmin.y), tmin.z);
    let dst_b = min(tmax.x, min(tmax.y, tmax.z)); 

    let dst_to_box = max(0.0, dst_a); 
    let dst_inside_box = max(0.0, dst_b - dst_to_box);
    return vec2(dst_to_box, dst_inside_box);
}

fn linearize_depth(d: f32, near: f32) -> f32 {
    return near / d;
}

const ITERATIONS: i32 = 500;

@fragment
fn fragment(
    mesh: MeshVertexOutput,
) -> @location(0) vec4<f32> {
    var depth = bevy_pbr::prepass_utils::prepass_depth(mesh.position, 0u);
    depth = linearize_depth(depth, 0.1);
    let v = normalize(mesh.world_position.xyz - view.world_position);

    let intersection = ray_box(view.world_position, v, settings.box_min, settings.box_max);
    let through_box = (min(intersection.x + intersection.y, depth) - intersection.x);

    let step_size = through_box / f32(ITERATIONS);

    
    var acc_pos = 0.0; 
    var acc_neg = 0.0;
    var cum_density = 0.0;
    for (var i = 0; i < ITERATIONS; i++) {
        let p = view.world_position + v * (intersection.x + f32(i) * step_size);
        let remapped_to_box = (p - settings.box_min) / (settings.box_max - settings.box_min);
         
        
        let density = textureSample(density_texture, density_sampler, remapped_to_box).x;
        let prob = density * density * settings.density_multiplier;

        let transmittance = exp(-cum_density);

        acc_pos += max(0.0, density) * transmittance; 
        acc_neg -= min(0.0, density) * transmittance;
        cum_density += prob * step_size;
    }

    let positive = acc_pos / (acc_pos + acc_neg + 0.1);
    let negative = acc_neg / (acc_pos + acc_neg + 0.1);

    let wave_color = mix(vec3(min(1.0, negative), 0.0, 0.0), vec3(0.0, 0.0, min(1.0, positive)), positive / max(1e-6, positive + negative));

    return vec4(wave_color, cum_density);
}