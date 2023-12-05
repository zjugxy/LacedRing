#pragma once
#include <unordered_set>
#include<unordered_map>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include<queue>
#include"Meshlet.h"
typedef OpenMesh::TriMesh_ArrayKernelT<>  MyMesh;


struct Meshlet_built
{
	std::unordered_set<uint> vertices;
	std::unordered_set<uint> faces;
	
};

enum SpecialTopo
{
	Normal,
	Surrouned,
	Include
};


struct EwireBuildSet {
	std::unordered_set<uint> borderset;
	std::unordered_set<uint> cornerset;
	std::unordered_set<uint> halfedgeset;
	std::vector<std::vector<uint>> ewires;//vertex wire
	std::vector<uint> ewiresid;
	std::vector<bool> ewirereverse;

	SpecialTopo sign = Normal;
};

class NewCluster
{
public:

	std::vector<Meshlet_built> mymeshlets;
	std::vector<std::vector<int>> oldmeshlets;
	//
	//std::vector<std::vector<uint>> Ewires;
	size_t nfaces;
	size_t nvertices;

	std::vector<std::vector<uint>> gloewires;
	std::vector<uint> lateaddress_surred;
	std::vector<uint> lateaddress_inc;

public:
	NewCluster(uint maxverts, uint maxtris, const MyMesh& mesh);

private:
	bool AddToMeshlet(int maxv, int maxt, Meshlet_built* curr, std::array<uint, 3> pnts,uint& faceid);
	void CalculateCenterRadius(const std::vector<vec3>& positions, float& radius, vec3& center);

	float ComputeScore(const Meshlet_built* curr, const vec3& center, const float radius,
		const vec3& normal, const std::array<vec3, 3>& tri, const std::array<uint, 3>& triidxs);
	bool IsMeshletFull(const Meshlet_built* curr,const int maxv,const int maxf);
	
	void PackOldmeshlets();
	void GenerateEwire(const MyMesh& mesh);

	uint CalculateVertexCnt(const MyMesh& mesh, uint vidx,const std::vector<uint>& meshletid);
	void BuildEwiresOfSingleMeshlet(EwireBuildSet& buildset, const MyMesh& mesh);


	void PushEdgeAtback(std::deque<uint>& edgewire, const MyMesh& mesh, std::unordered_set<uint>& edges, const std::unordered_set<uint> corners);
	void PushEdgeAtFront(std::deque<uint>& edgewire, const MyMesh& mesh, std::unordered_set<uint>& edges, const std::unordered_set<uint> corners);
	void PackAndCheck(std::vector<std::vector<uint>>& ewires , const std::unordered_set<uint>& corners, std::deque<uint>& edgewire, const MyMesh& mesh);
	void BuildGloEwires(std::vector<EwireBuildSet>& buildsets);
	void ReorderAndCheck(EwireBuildSet& buildset);
};

