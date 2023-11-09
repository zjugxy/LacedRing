#version 450 core
#extension GL_NV_mesh_shader : require

#define GROUP_SIZE 32

layout(local_size_x=GROUP_SIZE) in;
layout(max_vertices=96, max_primitives=32) out;
layout(triangles) out;


layout(std430, binding=1) readonly buffer vertex_data
{
    vec4 position[];
};


struct TS_meshlet
{
	uint vertex_begin;
	uint primitive_cnt;
	uint useless1;
	uint useless2;
	vec4 color;
};

layout(std430, binding=2) readonly buffer meshletdata
{
    TS_meshlet meshlets[];
};

uniform mat4 model;
uniform mat4 view;
uniform mat4 projection;

// mesh output
layout (location = 0) out myPerVertexData
{
  vec3 color;
} v_out[]; 



void main()
{
	int num=9;

		TS_meshlet meshdes = meshlets[gl_WorkGroupID.x];
		uint begin = meshdes.vertex_begin;
		uint irc = gl_LocalInvocationID.x*3;
		uint mi = gl_WorkGroupID.x;
		uint threadid = gl_LocalInvocationID.x;

	//每个线程处理一个图元，因为是triangle soup 处理三个顶点，添加三个index
	if(gl_LocalInvocationID.x<meshdes.primitive_cnt&&mi==0&&threadid==0){

		gl_MeshVerticesNV[gl_LocalInvocationID.x*3].gl_Position = projection*view*model*position[begin+irc];
		gl_MeshVerticesNV[gl_LocalInvocationID.x*3+1].gl_Position = projection*view*model*position[begin+1+irc];
		gl_MeshVerticesNV[gl_LocalInvocationID.x*3+2].gl_Position = projection*view*model*position[begin+2+irc];

		v_out[irc+0].color = vec3(meshdes.color);
		v_out[irc+1].color = vec3(meshdes.color);
		v_out[irc+2].color = vec3(meshdes.color); 
	
		gl_PrimitiveIndicesNV[gl_LocalInvocationID.x*3] = gl_LocalInvocationID.x*3;
		gl_PrimitiveIndicesNV[gl_LocalInvocationID.x*3+1] = gl_LocalInvocationID.x*3+1;
		gl_PrimitiveIndicesNV[gl_LocalInvocationID.x*3+2] = gl_LocalInvocationID.x*3+2;
	}


	if(gl_LocalInvocationID.x==0)
		gl_PrimitiveCountNV = meshdes.primitive_cnt;


}
   