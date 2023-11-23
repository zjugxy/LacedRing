#version 450 core
#extension GL_NV_mesh_shader : require
#define GROUP_SIZE 32

layout(local_size_x=GROUP_SIZE) in;
layout(max_vertices=64, max_primitives=126) out;//error may occur
layout(triangles) out;

layout(std430,binding =0)readonly buffer MeshletOff{
	unsigned int meshletoff[];
}

layout(std430,binding = 1)readonly buffer Geoinfo{
	vec4 position[];
}

layout(std430,binding = 2)readonly buffer MeshletDataPart1{
	unsigned char datapart1[];
}

layout(std430,binding = 3) readonly buffer MeshletDataPart2{
	unsigned int datapart2[];
}

layout(std430,binding =4)readonly buffer ExternalData{
	unsigned char externaldata[];
}


void main(){
	uint meshletid = gl_WorkGroupID.x;
	uint threadid = gl_LocalInvocationID.x;

	unsigned int offset1 = meshletoff[2*meshletid];
	unsigned int offset2 = meshletoff[2*meshletid+1];

	unsigned char




}