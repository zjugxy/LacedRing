#pragma once
#include<vector>
#include<unordered_set>
#include"NewCluster.h"

using uchar = unsigned char;
using uint = unsigned int;

struct MeshletBody {
	std::vector<uchar>left;
	std::vector<uchar> right;
	std::vector<uchar> irregular;
	std::vector<uint> vertex;
};

struct LWMeshlet
{
	MeshletBody InternalLW;
	std::vector<uint> ewireidx;
	std::vector<bool> reverse;
};

struct LWbuilder
{
	std::unordered_set<uint> vertices;
	std::unordered_set<uint> faces;

	std::unordered_set<uint> cornerset;
	std::unordered_set<uint> borderset;
};

struct MeshletDes
{
	uchar ewirenum;
	uchar color[3];
	uchar irrnum;
	uchar numvertex;
	uchar useless1;
	uchar useless2;
	uint ingeolocation;// 6+ 26
	uint inconlocation;
	std::vector<uint> exgeolocation;//1+5+26
	std::vector<uint> exconlocation;
};



class NewLWGenerator
{
public:
	std::vector<MeshletBody> gloEwires;
	std::vector<LWbuilder> internalbuilders;
	std::vector<LWMeshlet> targets;
	uint nummeshlet;
	std::vector<std::unordered_map<uint, uchar>> vidxmaps;

	//simple lw
	std::vector<Simple_meshlet> simplemeshlets;
	std::vector<vec4> geoinfo;
	std::vector<uchar> priminfo;
	// gpu lw
	std::vector<uint> DesLoc;
	std::vector<uint> Desinfo;

	std::vector<uint> newintercon;
	std::vector<uint> newextercon;
	std::vector<uchar> intercon;//inter: left,right,irregular
	std::vector<uchar> extercon;//exter: left,right
	std::vector<float> intergeo;
	std::vector<float> extergeo;



public:
	NewLWGenerator(const NewCluster& nclu);

	void FUNC(const MyMesh& mesh);

private:
	void InternalWireGenerator(const MyMesh& mesh);
	void LeftRightGenertor(const MyMesh& mesh);


private:
	void LevelPntGen(const std::unordered_set<uint>& bordset,const std::unordered_set<uint>& verset,
		std::vector<std::unordered_set<uint>>& lpnts,const MyMesh&mesh);
	void PushAtBack(std::vector<std::unordered_set<uint>>& lpnts, const MyMesh& mesh, std::deque<uint>& singlewire, std::unordered_set<uint>& pntstoremove);
	void PushAtFront(std::vector<std::unordered_set<uint>>& lpnts, const MyMesh& mesh, std::deque<uint>& singlewire, std::unordered_set<uint>& pntstoremove);
	void PackSimpleLaceWire(const MyMesh& mesh);

	void PackGPULW(const MyMesh& mesh);

	std::array<uchar, 4> DiscomposeUint(uint value);
	float PackChar4Float(uchar c0, uchar c1, uchar c2, uchar c3);
	uint PackChar4Uint(uchar c0, uchar c1, uchar c2, uchar c3);
};

