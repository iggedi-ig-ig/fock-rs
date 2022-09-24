#version 450

layout (binding = 0) uniform texture3D densityMap;
layout (binding = 1) uniform sampler smplr;

layout (location = 0) in vec2 fragCoord;
layout (location = 0) out vec4 fragColor;

layout (push_constant) uniform Uniforms {
    float dens_multiplier;
    float pos_x;
    float pos_y;
    float pos_z;
    float yaw;
    float pitch;
    float invAspect;
    uint nvox;
} uniforms;

mat2 rot2D(float a) {
    float s = sin(a);
    float c = cos(a);

    return mat2(c, -s, s, c);
}

void main() {
    vec3 ro = vec3(uniforms.pos_x, uniforms.pos_y, uniforms.pos_z);
    vec3 rd = normalize(vec3((fragCoord.xy - 0.5) * vec2(1.0, uniforms.invAspect), 0.5));

    rd.yz *= rot2D(uniforms.pitch);
    rd.xz *= rot2D(-uniforms.yaw);

    float accum_pos = 0.0;
    float accum_neg = 0.0;

    const float iters = 250.0;
    for (int i = 0; i < iters; i++) {
        vec3 pos = ro + rd * float(i) * 0.05;
        float val = texture(sampler3D(densityMap, smplr), pos * 0.25).r;

        accum_pos += max(0.0, val);
        accum_neg -= min(0.0, val);
    }

    fragColor = vec4(log(1.0 + vec3(accum_neg, 0.0, accum_pos) / iters * uniforms.dens_multiplier), 1.0);
}
