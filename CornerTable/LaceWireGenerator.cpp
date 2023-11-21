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


		//find start cycle vertex
		uint starthlidx;
		uint startidx = FindStartVertex(meshlet, mesh, boundvset,starthlidx);
		//
		//���ܻ���Ҫ����single internal pont�����
		//
		if (startidx == NOINTERNALVERTEX) {
			assert(boundvset.size() == meshlet.vertexnum);
			//����ӳ��
			//external idx --> meshlet limit idx map
			std::map<uint, short> gloidx2idx;
			short startidx = 0;
			for (auto& id : meshlet.externalwireids)
				for (int i = 0; i < Ewires[id].wire.size() - 1; i++)
					gloidx2idx[Ewires[id].wire[i]] = startidx++;



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
		else {
			//��startidx��ʼbuild loop
			//�����İ�
		}

	}

}

LaceWireGenerator::LaceWireGenerator()
{
}

uint LaceWireGenerator::FindStartVertex(const LaceWire_meshlet& meshlet, const MyMesh& mesh,const std::unordered_set<uint>& boundset,uint& hlidx){
	
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
					if (meshlet.faces.find(faceid) == meshlet.faces.end())std::cout << "error" << std::endl;
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
