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
    box_size: f32,
    pos_x: f32,
    pos_y: f32,
    pos_z: f32,
    yaw: f32,
    pitch: f32,
    invAspect: f32,
}

struct FragmentOutput {
    @location(0) fragColor: vec4<f32>,
}

@group(0) @binding(0) var densityMap: texture_3d<f32>;
@group(0) @binding(1) var smplr: sampler;
@group(0) @binding(2) var<storage> atomBuf: AtomBuffer;

var<push_constant> uniforms: Uniforms;

let iters: f32 = 350.0;

fn rayAABB(ro: vec3<f32>, rd: vec3<f32>, bmin: vec3<f32>, bmax: vec3<f32>) -> vec2<f32> {
    let inv_rd = 1.0 / rd;

    let t1 = (bmin - ro) * inv_rd;
    let t2 = (bmax - ro) * inv_rd;

    let tmin = max(max(min(t1.x, t2.x), min(t1.y, t2.y)), min(t1.z, t2.z));
    let tmax = min(min(max(t1.x, t2.x), max(t1.y, t2.y)), max(t1.z, t2.z));

    return vec2<f32>(tmin, tmax);
}

fn raySphere(ro: vec3<f32>, rd: vec3<f32>, center: vec3<f32>, radius: f32) -> vec2<f32> {
    let m = ro - center;
    let b = dot(m, rd);
    let c = dot(m, m) - radius * radius;

    if (c > 0.0 && b > 0.0f) { return vec2<f32>(0.0, -1.0); }

    let discr = b * b - c;
    if (discr < 0.0f) { return vec2<f32>(0.0, -1.0); }

    let r = sqrt(discr);
    return vec2<f32>(-b - r, -b + r);
}

fn rot2D(a: f32) -> mat2x2<f32> {
    let s = sin(a);
    let c = cos(a);

    return mat2x2<f32>(c, -s, s, c);
}

fn atom(ro: vec3<f32>, rd:vec3<f32>) -> vec4<f32> {
    var min_t = 1000.0;
    var color: vec3<f32> = vec3<f32>(0.0);
    for (var i = 0; i < atomBuf.size; i++) {
        let atom = atomBuf.atoms[i];
        let atom_pos = vec3<f32>(atom.pos_x, atom.pos_y, atom.pos_z);
        let atom_col = vec3<f32>(atom.red, atom.green, atom.blue);
        let atom_radius = atom.radius;
        let intersection = raySphere(ro, rd, atom_pos, atom_radius);

        if (intersection.x < intersection.y) {
            let t = intersection.x;

            if (t < min_t) {
                min_t = t;
                color = atom_col;
            }
        }
    }

    return vec4<f32>(color, min_t);
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

    // render molecule
    let result = atom(ro, rd);
    let col = result.rgb;
    let dist = result.w;

    // bounding volume intersection
    let box_size = uniforms.box_size;
    let intersection = rayAABB(ro, rd, vec3(-box_size * 0.5), vec3(box_size * 0.5));

    // ray through bounding volume
    var t_start = max(0.0, intersection.x);
    var t_end = min(dist, intersection.y);

    let ray_distance = t_end - t_start;
    let step_size = max(ray_distance, 1e-2) / iters;

    var acc_pos = 0.0;
    var acc_neg = 0.0;
    var acc_density = 0.0;
    for (var t = t_start; t < t_end && exp(-acc_density) > 1e-3; t += step_size) {
        let p = ro + rd * t;

        // wave function value at p
        let wave = textureSample(densityMap, smplr, p / box_size + 0.5).r;
        let prob = wave * wave * uniforms.dens_multiplier;

        // amount of reflected / absorbed intensity of the light at this position
        // probability of light reaching this far into the probability density

        let transmittance = exp(-acc_density);

        acc_pos += max(0.0, wave) * transmittance;
        acc_neg -= min(0.0, wave) * transmittance;

        acc_density += prob * step_size;
    }

    let positive = acc_pos / (acc_pos + acc_neg);
    let negative = acc_neg / (acc_pos + acc_neg);

    let waveCol = mix(vec3<f32>(min(1.0, negative), 0.0, 0.0), vec3<f32>(0.0, 0.0, min(1.0, positive)), positive / max(1e-6, positive + negative));
    let color = mix(col, waveCol, min(1.0, acc_density));

    return FragmentOutput(vec4<f32>(color, 1.0));
}
