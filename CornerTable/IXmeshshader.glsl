#version 450 core
#extension GL_NV_mesh_shader : require

#define GROUP_SIZE 32

layout(local_size_x=GROUP_SIZE) in;
layout(max_vertices=64, max_primitives=126) out;
layout(triangles) out;


layout(std430,binding = 1) readonly buffer geoinfo{
	vec4 position[];
};

layout(std430,binding = 2) readonly buffer verinfo{
	int veridx[];
};

layout(std430,binding = 3) readonly buffer priminfo{
	int primidx[];
};

struct IX_meshlet
{
	uint vertex_cnt;
	uint vertex_begin;
	uint primcnt;
	uint primbegin;
	vec4 color;
};

layout(std430,binding = 4) readonly buffer meshinfo{
	IX_meshlet meshlets[];
};

uniform mat4 model;
uniform mat4 view;
uniform mat4 projection;

// mesh output
layout (location = 0) out myPerVertexData
{
  vec3 color;
} v_out[]; //[max_vertices]


void main(){
	uint mi = gl_WorkGroupID.x;
	uint threadid = gl_LocalInvocationID.x;
	IX_meshlet meshlet = meshlets[mi];

	//每个线程处理顶点


		uint add = threadid;
		while(add<meshlet.vertex_cnt){
			int gloveridx = veridx[meshlet.vertex_begin+add];
			vec4 vergeo = position[gloveridx];
			gl_MeshVerticesNV[add].gl_Position = projection*view*model*vergeo;
			v_out[add].color = vec3(meshlet.color);
			add+=GROUP_SIZE;
		}

		add = threadid;
		while(add<meshlet.primcnt){
			uint start = meshlet.primbegin + add*3;
			gl_PrimitiveIndicesNV[add*3] = primidx[start];
			gl_PrimitiveIndicesNV[add*3+1] = primidx[start+1];
			gl_PrimitiveIndicesNV[add*3+2] = primidx[start+2];
			add+=GROUP_SIZE;
		}

		if(threadid==0)
			gl_PrimitiveCountNV = meshlet.primcnt;


}