#pragma once
#include"Meshlet.h"
#include<unordered_set>
#include<deque>

struct FinalMeshlet {
	unsigned char vertexcnt;
	unsigned char facecnt;
	unsigned char ewirecnt;
	unsigned char __meaningless;

	uint ptr_inter;
	std::vector<uint> ptrs_exter;
	//别说了，直接变成uint大军
};


class LaceWireGenerator
{
public:
	std::vector<ExternalWire> Ewires;
	std::vector<LaceWire_meshlet> meshlets;
	std::map<int, int> Dual2idx;
	std::vector<std::map<uint, short>> gloidx2idxvec;
	std::vector<std::map<short, uint>> idx2glo;

	std::vector<vec4> geoinfo;
	std::vector<unsigned char> priminfo;
	std::vector<Simple_meshlet> simplemeshlets;


	//GPU packer

	std::vector<uint> FinalMeshletData;
	std::vector<uint> meshletdata;

public:
	void InternalWireGeneraotr(const MyMesh& mesh);
	LaceWireGenerator();
	//分两步验证 internal wire is right?
	// meshlet triangle is right?


	//uint FindStartVertex(const LaceWire_meshlet& meshlet, const MyMesh& mesh,const std::unordered_set<uint>& boundset,uint& hlidx);
	short FindLeftOrRightIDX_VERSION0(uint start,uint end,const LaceWire_meshlet& meshlet, const MyMesh& mesh,  std::map<uint, short>& idxmap,short& leftorRight, std::set<uint>& facestoremove);
	short FindLeftOrRightIDX_VERSION1(uint start, uint end, const LaceWire_meshlet& meshlet, const MyMesh& mesh, std::map<uint, short>& idxmap, short& leftorRight,std::set<uint>&faces);
	void FindLeftAndRightIDX_VERSION2(uint start, uint end, const LaceWire_meshlet& meshlet, const MyMesh& mesh, std::map<uint, short>& idxmap, short& leftidx,short& Rightidx, std::set<uint>& faces);
	
	void LevelWireGenerator( LaceWire_meshlet& meshlet, const MyMesh& mesh, const std::unordered_set<uint>& boundset);
	void LevelOfontGeneratr(std::vector<std::set<uint>>& levelpnts,const LaceWire_meshlet& meshlet, const MyMesh& mesh, const std::unordered_set<uint>& boundset);
	uint GetstartPnts(std::vector<std::set<uint>>& levelpnts);
	void GenerateOneWire(std::deque<uint>& onewire ,int& remainnums,std::vector<std::set<uint>>& levelofPnts,const MyMesh& mesh);
	
	bool IsNextPnt(uint vidx, std::vector<std::set<uint>>& levelofPnts);

	void PackIntoGPUSimple(const MyMesh& mesh);
	void PackIntoGPU(const MyMesh& mesh);
	void SimpleCheckPrimIdx(int primbegin, const MyMesh& mesh, const LaceWire_meshlet& meshlet,int cntid);

	 uint Packunsignedchar(unsigned char& a, unsigned char& b, unsigned char& c, unsigned char d);
	// Interlacewire(std::vector<std::array<uint, 2>>& wirerecord, const MyMesh& mesh, LaceWire_meshlet& meshlet, const std::unordered_set<uint>& boundset);
	//bool Internextpntsearch(const MyMesh& mesh, const std::unordered_set<uint>& boundset, uint& startidx, uint& starthlidx, uint& nextidx, uint& nexthlidx,const std::set<uint>& pre);
};

