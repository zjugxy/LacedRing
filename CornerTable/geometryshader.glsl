#version 330 core

layout(triangles) in;
layout(triangle_strip,max_vertices=3) out;

//uniform int faceidtocolor[5002];
uniform isampler1D tex0;
uniform isampler1D tex1;

uniform bool UseCT;

in vec3 vFragPos[3];

out vec3 Normal;
out vec3 FragPos;
out vec3 basecolor;


void main(){
	int faceID = gl_PrimitiveIDIn;
	vec3 n = normalize(cross(vFragPos[1].xyz-vFragPos[0].xyz, vFragPos[2].xyz-vFragPos[0].xyz));
	
	int faceidtocolor = 0;
	if(UseCT)
		 faceidtocolor = texelFetch(tex0,faceID,0).r;
	else
		faceidtocolor = texelFetch(tex1,faceID,0).r;

	for(int i=0;i<3;i++){
		if(faceidtocolor==0)
			basecolor = vec3(0.6f,0.9f,1.0f);
		else if(faceidtocolor==1)
			basecolor = vec3(1.0f,0.82f,0.0f); 
		else if(faceidtocolor==2)
			basecolor = vec3(0.0f,0.47f,0.62f);
		else if(faceidtocolor==3)
			basecolor = vec3(0.0f,1.0f,1.0f);
		else if(faceidtocolor==4)
			basecolor = vec3(1.0f,0.0f,1.0f);
		else
			basecolor = vec3(1.0f,1.0f,1.0f);

		//basecolor = vec3(1.0f,0.0f,1.0f);
		gl_Position = gl_in[i].gl_Position;
		Normal = n;
		FragPos = vFragPos[i];
		EmitVertex();
	}

	EndPrimitive();
}