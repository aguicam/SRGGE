#version 330

smooth in vec3 N;
smooth in vec3 V;
smooth in vec2 vTexCoord;
out vec4 frag_color;
uniform int lod;
uniform sampler2D wall_tex;

void main (void) {

   frag_color =texture(wall_tex, vTexCoord);

}
