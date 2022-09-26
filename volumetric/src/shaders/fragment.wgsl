struct Atom {
    radius: f32,
    red: f32,
    green: f32,
    blue: f32,
    pos_x: f32,
    pos_y: f32,
    pos_z: f32,
}

struct AtomBuffer {
    size: i32,
    atoms: array<Atom>,
}

struct Uniforms {
    dens_multiplier: f32,
    pos_x: f32,
    pos_y: f32,
    pos_z: f32,
    yaw: f32,
    pitch: f32,
    invAspect: f32,
    nvox: u32,
}

struct FragmentOutput {
    @location(0) fragColor: vec4<f32>,  
}

struct MarchResult {
    dist: f32,
    color: vec3<f32>,
}

@group(0) @binding(0) var densityMap: texture_3d<f32>;
@group(0) @binding(1) var smplr: sampler;
@group(0) @binding(2) var<storage> atomBuf: AtomBuffer;

var<push_constant> uniforms: Uniforms;

let march_iters: i32 = 100;
let march_epsilon: f32 = 1e-2;

let iters: f32 = 200.0;
let step_size: f32 = 0.1;

fn rot2D(a: f32) -> mat2x2<f32> {
    let s = sin(a);
    let c = cos(a);

    return mat2x2<f32>(c, -s, s, c);
}

fn molecule(p: vec3<f32>) -> vec4<f32> {
    var closest: i32;
    var minDistSq: f32 = 1e12;
    for (var i = 0; i < atomBuf.size; i++) {
        let atom = atomBuf.atoms[i];
        let diff = vec3<f32>(atom.pos_x, atom.pos_y, atom.pos_z) - p;

        let distSq = dot(diff, diff);

        if (distSq < minDistSq) {
            minDistSq = distSq;
            closest = i;
        }
    }

    let atom = atomBuf.atoms[closest];
    let color = vec3<f32>(atom.red, atom.green, atom.blue);

    return vec4<f32>(color, sqrt(minDistSq) - atom.radius);
}

fn normal(p: vec3<f32>) -> vec3<f32> {
    let e = vec2<f32>(0.001, 0.0);
    var n: vec3<f32>;

    n.x = molecule(p + e.xyy).w - molecule(p - e.xyy).w;
    n.y = molecule(p + e.yxy).w - molecule(p - e.yxy).w;
    n.z = molecule(p + e.yyx).w - molecule(p - e.yyx).w;

    return normalize(n);
}

fn march(ro: vec3<f32>, rd: vec3<f32>) -> MarchResult {
    var dist: f32 = 0.0;

    for (var i = 0; i < march_iters; i++) {
        let p: vec3<f32> = ro + rd * dist;
        let res = molecule(p);

        dist += res.w;

        if (res.w < march_epsilon) {  return MarchResult(dist, res.rgb); }
        else if (dist > 50.0) { break; }
    }

    return MarchResult(dist, vec3<f32>(0.0));
}

@fragment
fn main(@location(0) fragCoord: vec2<f32>) -> FragmentOutput {
    let ro = vec3<f32>(uniforms.pos_x, uniforms.pos_y, uniforms.pos_z);
    var rd = normalize(vec3<f32>((fragCoord.xy - 0.5) * vec2<f32>(1.0, uniforms.invAspect), 0.5));

    // wgsl bad
    let rd_rotated_yz = rd.yz * rot2D(uniforms.pitch);
    rd.y = rd_rotated_yz.x;
    rd.z = rd_rotated_yz.y;

    let rd_rotated_xz = rd.xz * rot2D(-uniforms.yaw);
    rd.x = rd_rotated_xz.x;
    rd.z = rd_rotated_xz.y;

    let result = march(ro, rd);
    let dist = result.dist;
    var col = vec3<f32>(result.color);

    if (dist < 50.0) {
        let n = normal(ro + rd * dist);
        let diff = dot(-n, rd);

        col *= diff * 0.9 + 0.1;
    }

    var pos: f32 = 0.0;
    var neg: f32 = 0.0;

    var was_in = false;
    for (var i = 0; i < i32(iters); i++) {
        let d = f32(i) * step_size;

        let p = ro + rd * min(d, dist + march_epsilon);

        let inside = abs(p.x) < 10.0 && abs(p.y) < 10.0 && abs(p.z) < 10.0;
        if (was_in && !inside) {
            break;
        } else if (!was_in && inside) {
            was_in = true;
        }

        let val = textureSample(densityMap, smplr, p / 20.0 + 0.5).r;

        pos += max(val, 0.0);
        neg -= min(val, 0.0);

        if (d > dist + march_epsilon) { break; }
    }

    let wave = uniforms.dens_multiplier * (neg + pos) / iters;
    let prob = wave * wave;

    let color = mix(col, vec3<f32>(neg, 0.0, pos), prob);
    return FragmentOutput(vec4<f32>(color, 1.0));
}
