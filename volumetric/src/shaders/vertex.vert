#version 430

layout (location = 0) out vec2 fragCoord;

void main() {
    // 0---1
    // | / |
    // 2---3

    vec2 soos = vec2(float(gl_VertexIndex & 1), float(gl_VertexIndex >> 1));
    fragCoord = vec2(soos.x, soos.y);
    gl_Position = vec4(soos * 2 - 1, 0.0, 1.0);
}