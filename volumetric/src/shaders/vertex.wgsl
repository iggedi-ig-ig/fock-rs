struct VertexOutput {
    @location(0) fragCoord: vec2<f32>,
    @builtin(position) member: vec4<f32>,
}

var<private> fragCoord: vec2<f32>;
var<private> gl_VertexIndex: u32;
var<private> gl_Position: vec4<f32>;

fn main_1() {
    var soos: vec2<f32>;

    let _e2 = gl_VertexIndex;
    let _e7 = gl_VertexIndex;
    soos = vec2<f32>(f32((_e2 & u32(1))), f32((_e7 >> u32(1))));
    let _e14 = soos;
    let _e16 = soos;
    fragCoord = vec2<f32>(_e14.x, _e16.y);
    let _e20 = soos;
    let _e27 = ((_e20 * f32(2)) - vec2<f32>(f32(1)));
    gl_Position = vec4<f32>(_e27.x, _e27.y, 0.0, 1.0);
    return;
}

@vertex 
fn main(@builtin(vertex_index) param: u32) -> VertexOutput {
    gl_VertexIndex = param;
    main_1();
    let _e5 = fragCoord;
    let _e7 = gl_Position;
    return VertexOutput(_e5, _e7);
}
