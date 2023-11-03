#version 330 core

layout(triangles) in;
layout(triangle_strip,max_vertices=3) out;

//uniform int faceidtocolor[5002];
uniform isampler1D tex0;
uniform isampler1D tex1;
uniform isampler1D tex2;

uniform int n_face;
uniform bool UseCT;

in vec3 vFragPos[3];

out vec3 Normal;
out vec3 FragPos;
out vec3 basecolor;


void main(){
	int faceID = gl_PrimitiveIDIn;
	vec3 n = normalize(cross(vFragPos[1].xyz-vFragPos[0].xyz, vFragPos[2].xyz-vFragPos[0].xyz));

	float texCoord = float(faceID) / float(n_face - 1);
    vec3 textureColor = texture(tex2, texCoord).rgb;

	for(int i=0;i<3;i++){
		gl_Position = gl_in[i].gl_Position;
		Normal = n;
		FragPos = vFragPos[i];
		EmitVertex();
	}

	EndPrimitive();
}