#version 450 core
#extension GL_NV_mesh_shader : require
#define GROUP_SIZE 32

#define WIREEND 0xFF
#define EMPTYWIRE 0xFFFFFFFF

layout(local_size_x=GROUP_SIZE) in;
layout(max_vertices=64, max_primitives=126) out;//error may occur
layout(triangles) out;

layout(std430,binding = 1)readonly buffer PositionInfo{
	uint position[];
}

layout(std430,binding = 2)readonly buffer Desinfo{
	uint Desinfo[];
}

layout(std430,binding=3)readonly buffer Intermeshdata{
	float intermesh[];
}

layout(std430,binding=4)readonly buffer Extermeshdata{
	float extermesh[];
}


void main(){

	// unpack intermeshlet
	uint mi = gl_WorkGroupID.x;
	uint threadid = gl_LocalInvocationID.x;
	uint interlocation = position[mi];

	uint* Desinfostart = &Desinfo[interlocation];
	uint intermeshposition = Desinfostart[1];

	float* intermeshstart = &intermesh[intermeshposition];

	uint vertexcnt = (floatBitsToUint(intermeshstart[0])>>24)&0xFF;
	uint colorvertexoffset = (floatBitsToUint(intermeshstart[0])>>16)&0xFF;
	uint ircnt = (floatBitsToUint(intermeshstart[0])>>8)&0xFF;
	vec3 color = vec3(intermeshstart[colorvertexoffset],
					intermeshstart[colorvertexoffset+1],
					intermeshstart[colorvertexoffset+2])
	
	uint intprimcnt = vertexcnt*2-2+ircnt;

	uint add = threadid;
	while(add<intprimcnt){
		//left
		if(add<vertexcnt-1){
			gl_PrimitiveIndicesNV[add*3] = add%(vertexcnt);
			gl_PrimitiveIndicesNV[add*3+1] = (add+1)%(vertexcnt);
			gl_PrimitiveIndicesNV[add*3] = (intermeshstart[ (3+add%(vertexcnt))/4 ]>>(8*(3+add%(vertexcnt))%4))&0xFF;
		}
		//right
		else if(add<vertexcnt*2 -2){
			gl_PrimitiveIndicesNV[add*3] = add%(vertexcnt);
			gl_PrimitiveIndicesNV[add*3+1] = (add+1)%(vertexcnt);
			gl_PrimitiveIndicesNV[add*3] = (intermeshstart[ (3+add%(vertexcnt))/4 ]>>(8*(3+add%(vertexcnt))%4))&0xFF;
		}
		//irr
	
	}




}