#version 330

smooth in vec2 vTexCoord;
out vec4 frag_color;

uniform sampler2D floor_tex;

void main (void) {

   frag_color =texture(floor_tex, vTexCoord);

}
