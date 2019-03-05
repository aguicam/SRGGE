#version 330

layout (location = 0) in vec3 vert;
layout (location = 1) in vec3 normal;
layout (location = 2) in vec2 TexCoord;

uniform mat4 u_projection;
uniform mat4 u_view;
uniform mat4 u_model;
uniform mat3 u_normal_matrix;


smooth out vec3 N;
smooth out vec3 V;
smooth out vec2 vTexCoord;

void main(void)  {
    vec4 view_vertex = u_view * u_model * vec4(vert, 1);
    V = view_vertex.xyz;
    N = normalize(u_normal_matrix * normal);
    vTexCoord=TexCoord;
    gl_Position = u_projection * view_vertex;
}

