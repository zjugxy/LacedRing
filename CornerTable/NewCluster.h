#pragma once
#include <unordered_set>
#include<unordered_map>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include<queue>
#include<Eigen/Dense>
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

enum TopoHeal {
	Unchecked,
	Correct,
	Torus,
	SharpConnect
};

struct ShapeHeal {
	std::unordered_set<uint> halfedgeset;
	TopoHeal topotype = Unchecked;
	std::vector<uint> specialpnts;
	std::vector<uint> newfacestart;
};

class NewCluster
{
public:

	std::vector<Meshlet_built> mymeshlets;
	std::vector<std::vector<uint>> oldmeshlets;
	std::vector<ShapeHeal> detecters;
	//
	//std::vector<std::vector<uint>> Ewires;
	size_t nfaces;
	size_t nvertices;

	size_t maxf;
	size_t maxv;

	std::vector<std::vector<uint>> gloewires;
	std::vector<uint> lateaddress_surred;
	std::vector<uint> lateaddress_inc;
	std::vector<EwireBuildSet> buildsets;

	//
	std::vector<vec3> lines;


public:
	NewCluster(uint maxverts, uint maxtris, const MyMesh& mesh,Meshlets meshlets,const std::string& filename);

private:
	bool AddToMeshlet(int maxv, int maxt, Meshlet_built* curr, std::array<uint, 3> pnts,uint& faceid);
	void CalculateCenterRadius(const std::vector<vec3>& positions, float& radius, vec3& center);

	float ComputeScore(const Meshlet_built* curr, const vec3& center, const float radius,
		const vec3& normal, const std::array<vec3, 3>& tri, const std::array<uint, 3>& triidxs);

	bool IsMeshletFull(const Meshlet_built* curr,const int maxv,const int maxf);
	
	void PackOldmeshlets();
	void GenerateEwire(const MyMesh& mesh);

	void DeleteSmallMesh(const MyMesh& mesh);

	void FirstDeleteMesh(const MyMesh& mesh);

	void SimulateAnneal(const MyMesh& mesh,float threshold,float decreaserate,int cntloop);

	uint CalculateVertexCnt(const MyMesh& mesh, uint vidx,const std::vector<uint>& meshletid);

	void BuildEwiresOfSingleMeshlet(EwireBuildSet& buildset, const MyMesh& mesh);


	void PushEdgeAtback(std::deque<uint>& edgewire, const MyMesh& mesh, std::unordered_set<uint>& edges, const std::unordered_set<uint> corners);
	void PushEdgeAtFront(std::deque<uint>& edgewire, const MyMesh& mesh, std::unordered_set<uint>& edges, const std::unordered_set<uint> corners);
	void PackAndCheck(std::vector<std::vector<uint>>& ewires , const std::unordered_set<uint>& corners, std::deque<uint>& edgewire, const MyMesh& mesh);
	void BuildGloEwires(std::vector<EwireBuildSet>& buildsets);
	void ReorderAndCheck(EwireBuildSet& buildset);

	bool SharpTriCut(const MyMesh& mesh);
	void fitPlane(const Eigen::MatrixXd& points, Eigen::Vector3d& centroid, Eigen::Vector3d& normal);


	void TryShapeHeal(const MyMesh& mesh);
	bool ShapeErrorDetect(const MyMesh& mesh, uint mid);

	void TorusHeal(const MyMesh& mesh, uint mid);
	void SharpConnectHeal(const MyMesh& mesh, uint mid);
	void NewSharpConnectHeal(const MyMesh& mesh, uint mid);

	void GrowFaceSet(const MyMesh& mesh, std::unordered_set<uint>& newfaceset,
		std::unordered_set<uint>& faceset, uint seedface, uint cutlimit);
	void GrowFaceSet(const MyMesh& mesh, std::unordered_set<uint>& newfaceset,
		std::unordered_set<uint>& faceset, std::vector<uint> seedfaces, uint cutlimit);

	void GenerateTorusNewFaceStart(const MyMesh& mesh, const std::vector<std::vector<uint>>& loops, uint mid);
	
	std::vector<uint> BFS(const MyMesh& mesh, uint startidx,const std::unordered_set<uint>& targets,
		const std::unordered_set<uint>& faceset);

	double point2PlaneDist(const Eigen::Vector3d& point, const Eigen::Vector3d& centroid, const Eigen::Vector3d& normal);
	double EvaluateCost(const Eigen::Vector3d& centroid, const Eigen::Vector3d& normal, const std::vector<uint>& meshletid, uint faceid, const MyMesh& mesh);

	std::vector<uint> IsNodeBound(const std::vector<uint>& facetomeshlet, const MyMesh& mesh, int triid,int nodeid,uint& vertoerase);

	bool SimpleInsert(const std::vector<uint>& boundnodes, uint triid, const MyMesh& mesh,
		std::vector<uint>& facetomeshlet,uint rawnode,uint vertoerase);


	bool ComplexInsert(uint nodeinserted, uint triid, const MyMesh& mesh,
		std::vector<uint>& facetomeshlet, uint rawnode, uint vertoerase,std::unordered_set<uint> nodeoutsets);

	uint getVerCommon(const MyMesh& mesh, std::vector<uint> tris);

	void TestMeshlet(uint meshletid,const MyMesh& mesh);

	void PackBoundTris(const MyMesh& mesh, const Meshlet_built& meshlet, std::vector<uint>& tristopack, uint mid,const std::vector<uint>& facetomeshlet);

	void CalculateNormal(std::map<uint, vec3> normals, vec3& normal);


	void EvaluateMeshlet(const MyMesh& mesh);

	float FitsquareEnergy(const std::vector<Eigen::Vector3d>& points, Eigen::Vector3d& normal, Eigen::Vector3d& center);

	float FitsquareEnergy(const std::map<uint, Eigen::Vector3d>& positions, Eigen::Vector3d& normal, Eigen::Vector3d& center);

	float calculateTrianglearea(float a, float b, float c);

	float ComputeSimulateEnergy(float area, float length,const std::map<uint, Eigen::Vector3d>& positions);
};

