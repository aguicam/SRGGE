#version 330

smooth in vec3 N;
out vec4 frag_color;

void main (void) {

   frag_color =vec4(1,1,1,1)*N.z;

}
