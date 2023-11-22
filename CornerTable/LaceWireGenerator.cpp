#include "LaceWireGenerator.h"
#include<set>
#include<unordered_set>
#define NOINTERNALVERTEX 0xFFFFFFFF
#define NOTINI 0
#define RIGHT 1
#define LEFT 2


void LaceWireGenerator::InternalWireGeneraotr(const MyMesh& mesh)
{
	for (auto& meshlet : meshlets) {
		//���ɵ���meshlet�ڲ���internal wire
		
		std::unordered_set<uint> boundvset;
		for (auto& id : meshlet.externalwireids)
			for (auto& elem : Ewires[id].wire)
				boundvset.insert(elem);

		std::vector<uint> interidxseq;
		//
		meshlet.externalLR.resize(meshlet.externalwireids.size());
		for (auto& elem : meshlet.externalLR)elem = NOTINI;
		
		std::vector<short>& lrref = meshlet.externalLR;

		//********************
		//��Ҫ����һ��face���¼ T2 T1 T0������
		//********************

		if (boundvset.size() == meshlet.vertexnum) {
			//����ӳ��
			//external idx --> meshlet limit idx map
			std::map<uint, short> gloidx2idx;
			short mapidx = 0;
			for (auto& id : meshlet.externalwireids)
				for (int i = 0; i < Ewires[id].wire.size() - 1; i++)
					gloidx2idx[Ewires[id].wire[i]] = mapidx++;



			//����һ��wireֻ��Ҫ�洢 wire.size()-1 �� left right corner
			for (auto& id : meshlet.externalwireids)
				for (int i = 0; i < Ewires[id].wire.size() - 1; i++) {
					uint startidx = Ewires[id].wire[i];
					uint endidx = Ewires[id].wire[i+1];
					short LorR = NOTINI;
					short idx = FindLeftOrRightIDX_VERSION0(startidx, endidx, meshlet, mesh, gloidx2idx,LorR);
					if (lrref[id] == NOTINI)lrref[id] = LorR;
					else if (lrref[id] != LorR)std::cout << "error in find left or Right" << std::endl;

					if (LorR == RIGHT)Ewires[id].right.push_back(idx);
					else Ewires[id].left.push_back(idx);
				}

		}
		else if((boundvset.size() + 1)==meshlet.vertexnum){
			//��startidx��ʼbuild loop
			//�����İ�
			for(auto& ver:meshlet.vertexs)
				if (boundvset.find(ver) == boundvset.end()) {
					meshlet.interwire.wire.push_back(ver);
					break;
				}
			//mapӳ��
			std::map<uint, short> gloidx2idx;
			short mapidx = 0;
			gloidx2idx[meshlet.interwire.wire[0]] = mapidx++;
			for (auto& id : meshlet.externalwireids)
				for (int i = 0; i < Ewires[id].wire.size() - 1; i++)
					gloidx2idx[Ewires[id].wire[i]] = mapidx++;

			for (auto& id : meshlet.externalwireids)
				for (int i = 0; i < Ewires[id].wire.size() - 1; i++) {
					uint startidx = Ewires[id].wire[i];
					uint endidx = Ewires[id].wire[i + 1];
					short LorR = NOTINI;
					short idx = FindLeftOrRightIDX_VERSION0(startidx, endidx, meshlet, mesh, gloidx2idx, LorR);
					if (lrref[id] == NOTINI)lrref[id] = LorR;
					else if (lrref[id] != LorR)std::cout << "error in find left or Right" << std::endl;

					if (LorR == RIGHT)Ewires[id].right.push_back(idx);
					else Ewires[id].left.push_back(idx);
				}

		}
		else {
			//find start cycle vertex

			//first elem vertexidx second elem halfedge idx
			std::vector<std::array<uint, 2>> interwire;
			//�����������ɵĲ�����һ����������һ����
			Interlacewire(interwire, mesh, meshlet, boundvset);

		}

	}

}

LaceWireGenerator::LaceWireGenerator()
{
}

uint LaceWireGenerator::FindStartVertex(const LaceWire_meshlet& meshlet, const MyMesh& mesh,const std::unordered_set<uint>& boundset,uint& hlidx){
	//�ڶ���ʱ���ִ���
	for (auto& id : meshlet.externalwireids)
		for (int i = 0; i < Ewires[i].wire.size() - 1; i++) {
			int startvid = Ewires[i].wire[i];
			int endvid = Ewires[i].wire[i+1];

			OpenMesh::VertexHandle vh0(startvid);

			OpenMesh::VertexHandle vh1(endvid);

			for (auto it = mesh.cvoh_begin(vh0); it != mesh.cvoh_end(vh0); ++it) 
				//�ҵ�һ���߽�out half edge
				if (it->to().idx() == endvid) {
					int faceid = it->face().idx();
					
					//������start -->  end���ұ�
					if (meshlet.faces.find(faceid) != meshlet.faces.end()) {
						uint vidx = it->next().to().idx();
						if (boundset.find(vidx) == boundset.end()) {
							hlidx = it->idx();
							return vidx;
						}
						else
							continue;
					}
					//�����
					faceid = it->opp().face().idx();
					if (meshlet.faces.find(faceid) == meshlet.faces.end())
						std::cout << "error" << std::endl;
					uint vidx = it->opp().next().to().idx();
					if (boundset.find(vidx) == boundset.end()) {
						hlidx = it->opp().idx();
						return vidx;
					}
					else
						continue;
				}

		}
	return NOINTERNALVERTEX;

}

//���versionֻ������û��internal pnt�����
short LaceWireGenerator::FindLeftOrRightIDX_VERSION0(uint start, uint end, const LaceWire_meshlet& meshlet, const MyMesh& mesh, std::map<uint, short>& idxmap, short& leftorRight)
{
	OpenMesh::VertexHandle vh0(start);
	OpenMesh::VertexHandle vh1(end);

	for (auto it = mesh.cvoh_begin(vh0); it != mesh.cvoh_end(vh0); ++it)
		//�ҵ�һ���߽�out half edge
		if (it->to().idx() == end) {
			int faceid = it->face().idx();
			//�������ұ�
			if (meshlet.faces.find(faceid) != meshlet.faces.end()) {
				uint idx = it->next().to().idx();
				leftorRight = RIGHT;
				return idxmap[idx];
			}
			else {
			//face��wire���
				faceid = it->opp().face().idx();
				if (meshlet.faces.find(faceid) == meshlet.faces.end())std::cout << "error" << std::endl;
				uint idx = it->opp().next().to().idx();
				leftorRight = LEFT;
				return idxmap[idx];
			}
		}
}

void LaceWireGenerator::Interlacewire(std::vector<std::array<uint, 2>>& wirerecord, const MyMesh& mesh, LaceWire_meshlet& meshlet, const std::unordered_set<uint>& boundset)
{
	uint starthlidx;
	uint startidx = FindStartVertex(meshlet, mesh, boundset, starthlidx);
	uint nextidx, nexthlidx;
	uint intervertexnum = meshlet.vertexnum - boundset.size();
	wirerecord.push_back(std::array<uint, 2>{ startidx, 0});
	meshlet.interwire.wire.push_back(startidx);
	std::set<uint> alreadyin{ startidx };

	while (Internextpntsearch(mesh,boundset,startidx,starthlidx,nextidx,nexthlidx,alreadyin)) {
		if (meshlet.vertexs.find(nextidx) == meshlet.vertexs.end())
			std::cout << "error in occur outside pnts";

		auto edgehl = mesh.halfedge_handle(nexthlidx);
		edgehl = mesh.prev_halfedge_handle(edgehl);
		edgehl = mesh.opposite_halfedge_handle(edgehl);
		if (mesh.to_vertex_handle(edgehl).idx() != nextidx)
			std::cout << "error" << std::endl;
		if(mesh.from_vertex_handle(edgehl).idx()!=startidx)
			std::cout << "error" << std::endl;
		std::array<uint, 2> temp;
		temp[0] = nextidx; temp[1] = edgehl.idx();


		wirerecord.push_back(temp);
		alreadyin.insert(nextidx);
		meshlet.interwire.wire.push_back(nextidx);

		startidx = nextidx;
		starthlidx = nexthlidx;

	}

	if ((alreadyin.size() + boundset.size()) != meshlet.vertexnum) {
		std::cout << "not a full ring" << std::endl;
	}

}

bool LaceWireGenerator::Internextpntsearch(const MyMesh& mesh, const std::unordered_set<uint>& boundset, uint& startidx, uint& starthlidx, uint& nextidx, uint& nexthlidx, const std::set<uint>& pre)
{
	auto hl = mesh.halfedge_handle(starthlidx);
	
	do {
		hl = mesh.next_halfedge_handle(hl);
		hl = mesh.opposite_halfedge_handle(hl);
		hl = mesh.next_halfedge_handle(hl);
		auto vh = mesh.to_vertex_handle(hl);
		if (	(boundset.find(vh.idx()) == boundset.end())
				&&//vh.idx�Ȳ���boundset�ֲ���preset�м�
			(pre.find(vh.idx())==pre.end())	) {
			nextidx = vh.idx();
			nexthlidx = mesh.prev_halfedge_handle(hl).idx();
			return true;
		}


	} while (hl.idx() != starthlidx);
		return false;
}
