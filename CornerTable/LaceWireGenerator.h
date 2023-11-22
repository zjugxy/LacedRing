#pragma once
#include"Meshlet.h"
#include<unordered_set>
#include<deque>
class LaceWireGenerator
{
public:
	std::vector<ExternalWire> Ewires;
	std::vector<LaceWire_meshlet> meshlets;
	std::map<int, int> Dual2idx;

public:
	void InternalWireGeneraotr(const MyMesh& mesh);
	LaceWireGenerator();
	//分两步验证 internal wire is right?
	// meshlet triangle is right?


	//uint FindStartVertex(const LaceWire_meshlet& meshlet, const MyMesh& mesh,const std::unordered_set<uint>& boundset,uint& hlidx);
	short FindLeftOrRightIDX_VERSION0(uint start,uint end,const LaceWire_meshlet& meshlet, const MyMesh& mesh,  std::map<uint, short>& idxmap,short& leftorRight);
	void LevelWireGenerator( LaceWire_meshlet& meshlet, const MyMesh& mesh, const std::unordered_set<uint>& boundset);
	void LevelOfontGeneratr(std::vector<std::set<uint>>& levelpnts,const LaceWire_meshlet& meshlet, const MyMesh& mesh, const std::unordered_set<uint>& boundset);
	uint GetstartPnts(std::vector<std::set<uint>>& levelpnts);
	void GenerateOneWire(std::deque<uint>& onewire ,int& remainnums,std::vector<std::set<uint>>& levelofPnts,const MyMesh& mesh);
	
	bool IsNextPnt(uint vidx, std::vector<std::set<uint>>& levelofPnts);
	// Interlacewire(std::vector<std::array<uint, 2>>& wirerecord, const MyMesh& mesh, LaceWire_meshlet& meshlet, const std::unordered_set<uint>& boundset);
	//bool Internextpntsearch(const MyMesh& mesh, const std::unordered_set<uint>& boundset, uint& startidx, uint& starthlidx, uint& nextidx, uint& nexthlidx,const std::set<uint>& pre);
};

