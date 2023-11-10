#pragma once
#include"Viewer/viewermath.h"
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
typedef OpenMesh::TriMesh_ArrayKernelT<>  MyMesh;
typedef std::vector<std::vector<int>> Meshlets;
using uint = unsigned int;
//一个SSBO，其中装载了所有三角形序列

struct TS_meshlet
{
	uint vertex_begin;
	uint primitive_cnt;//cnt *3
	uint useless1;
	uint useless2;
	//range from    [vertex_begin,vertex_begin+primitive_cnt*3)
	vec4 color; // add one useless color float to fit std140
};

void TS_MeshletLoad(std::vector<TS_meshlet>& loader, MyMesh mesh, Meshlets meshlets,std::vector<vec4>& geoinfo);


struct IX_meshlet
{
	uint vertex_cnt;
	uint vertex_begin;
	uint primcnt;
	uint primbegin;
	vec4 color;
};

void IX_meshletLoad(std::vector<IX_meshlet>& loader, MyMesh mesh, Meshlets meshlets,std::vector<vec4>& geoinfo,std::vector<int>& vertexidx, std::vector<int>& primidx);

struct SC_meshlet {
	uint vertex_cnt;
	uint vertex_begin;
	uint primcnt;
	uint primbegin;
	vec4 color;

};

void SC_meshletLoad(std::vector<SC_meshlet>& loader, MyMesh mesh, Meshlets meshlets, std::vector<vec4>& geoinfo, std::vector<int>& primidx);

struct LW_meshlet
{
	uint internalwireid;
	std::vector<uint> externalwireid;
};
//
struct LWIDarray {
	std::vector<uint> prefixsum;
};

struct LacedWrie
{
	std::vector<uint> left;
	std::vector<uint> right;
	std::vector<vec3> geoinfo;
};

//需要一个前缀和数组标记每一个meshlet的start