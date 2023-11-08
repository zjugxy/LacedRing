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
	//range from    [vertex_begin,vertex_begin+primitive_cnt*3)
	vec3 color;
};

void MeshletLoad(std::vector<TS_meshlet>& loader, MyMesh mesh, Meshlets meshlets,std::vector<vec3>& geoinfo);

