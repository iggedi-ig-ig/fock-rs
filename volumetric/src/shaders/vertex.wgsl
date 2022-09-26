struct VertexOutput {
    @location(0) fragCoord: vec2<f32>,
    @builtin(position) position: vec4<f32>,
}

@vertex 
fn main(@builtin(vertex_index) vertex_index: u32) -> VertexOutput {
    let uv = vec2<f32>(f32(vertex_index & 1u), f32(vertex_index >> 1u));
    return VertexOutput(uv, vec4(2.0 * uv - 1.0, 0.0, 1.0));
}
