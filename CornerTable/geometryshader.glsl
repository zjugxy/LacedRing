#version 330 core

layout(triangles) in;
layout(triangle_strip,max_vertices=3) out;

uniform int faceidtocolor[320];

in vec3 vFragPos[3];

out vec3 Normal;
out vec3 FragPos;
out vec3 basecolor;


void main(){
	int faceID = gl_PrimitiveIDIn;
	vec3 n = normalize(cross(vFragPos[1].xyz-vFragPos[0].xyz, vFragPos[2].xyz-vFragPos[0].xyz));

	for(int i=0;i<3;i++){
		if(faceidtocolor[faceID]==0)
			basecolor = vec3(1.0f,1.0f,1.0f);//黄色
		else if(faceidtocolor[faceID]==1)
			basecolor = vec3(1.0f,0.0f,0.0f);  //对应T1为红色
		else if(faceidtocolor[faceID]==2)
			basecolor = vec3(0.0f,0.0f,1.0f);
		else if(faceidtocolor[faceID]==3)
			basecolor = vec3(0.0f,1.0f,1.0f);
		else if(faceidtocolor[faceID]==4)
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