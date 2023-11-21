#pragma once
#include"Meshlet.h"
#include<unordered_set>
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


	uint FindStartVertex(const LaceWire_meshlet& meshlet, const MyMesh& mesh,const std::unordered_set<uint>& boundset,uint& hlidx);
	short FindLeftOrRightIDX_VERSION0(uint start,uint end,const LaceWire_meshlet& meshlet, const MyMesh& mesh,  std::map<uint, short>& idxmap,short& leftorRight);
};

