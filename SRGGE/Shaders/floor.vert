#version 330

layout (location = 0) in vec3 vert;
layout (location = 1) in vec3 normal;
layout (location = 2) in vec2 TexCoord;

uniform mat4 u_projection;
uniform mat4 u_view;
uniform mat4 u_model;


smooth out vec2 vTexCoord;

void main(void)  {
    vec4 view_vertex = u_view * u_model * vec4(vert, 1);
    vTexCoord=TexCoord;
    gl_Position = u_projection * view_vertex;
}

