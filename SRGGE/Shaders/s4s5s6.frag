#version 330

smooth in vec3 N;
smooth in vec3 V;
out vec4 frag_color;
flat in vec3 light_pos;
uniform int lod;

uniform bool colored;
void main (void) {

    vec4 colorF=vec4(1,1,1,1)*N.z;
if(colored){
   if(lod==0){
       colorF=vec4(1,0,0,1)*N.z;
   }
   if(lod==1){
       colorF=vec4(0,1,0,1)*N.z;
   }
   if(lod==2){
       colorF=vec4(0,0,1,1)*N.z;
   }
   if(lod==3){
       colorF=vec4(1,0.5,0.12,1)*N.z;
   }
   if(lod==4){
       colorF=vec4(1,0,1,1)*N.z;
   }
}
   frag_color =colorF;

}
