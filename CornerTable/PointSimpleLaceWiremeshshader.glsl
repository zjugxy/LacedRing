#version 450 core
#extension GL_NV_mesh_shader : require

#define GROUP_SIZE 32

layout(local_size_x=GROUP_SIZE) in;
layout(max_vertices=64, max_primitives=150) out;
layout(points) out;


layout(std430,binding = 1) readonly buffer geoinfo{
	vec4 position[];
};

layout(std430,binding = 2) readonly buffer priminfo{
	unsigned char primidx[];
};

struct Simple_meshlet
{
	uint vertex_cnt;
	uint vertex_begin;
	uint primcnt;
	uint primbegin;
	vec4 color;
};

layout(std430,binding = 3) readonly buffer meshinfo{
	Simple_meshlet meshlets[];
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
	Simple_meshlet meshlet = meshlets[mi];
	if(mi==1){
	//ÿ���̴߳�����


		uint add = threadid;
		while(add<meshlet.vertex_cnt){
			uint gloveridx = meshlet.vertex_begin+add;
			vec4 vergeo = position[gloveridx];
			gl_MeshVerticesNV[add].gl_Position = projection*view*model*vergeo;
			v_out[add].color = vec3(meshlet.color);
			add+=GROUP_SIZE;
		}

		for(int i=0;i+threadid<meshlet.vertex_cnt;i+=GROUP_SIZE){
			gl_PrimitiveIndicesNV[i+threadid] = i+threadid;
		}

		if(threadid==0)
			gl_PrimitiveCountNV = meshlet.vertex_cnt;
	}
}