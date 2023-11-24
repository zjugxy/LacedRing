#include "LaceWireGenerator.h"
#include<set>
#include<unordered_set>
#include<deque>
#include<iostream>
#include<random>

#define NOINTERNALVERTEX 0xFFFFFFFF
#define NEXTWIRE 0xFFFFFFFF
#define EXTERNALWIREEND 0xFFFF

#define NOTINI 0
#define RIGHT 1
#define LEFT 2


void LaceWireGenerator::InternalWireGeneraotr(const MyMesh& mesh)
{
	int cntid = -1;
	gloidx2idxvec.resize(meshlets.size());
	for (auto& meshlet : meshlets) {
		//生成单个meshlet内部的internal wire
		cntid++;

		if (cntid == 72)
			std::cout << "the external wireid is wrong" << std::endl;

		std::unordered_set<uint> boundvset;
		for (auto& id : meshlet.externalwireids)
			for (auto& elem : Ewires[id].wire)
				boundvset.insert(elem);
		//
		meshlet.externalLR.resize(meshlet.externalwireids.size());
		for (auto& elem : meshlet.externalLR)elem = NOTINI;

		std::vector<short>& lrref = meshlet.externalLR;

		//********************
		//需要建立一个face表记录 T2 T1 T0并更新
		//********************

		if (boundvset.size() == meshlet.vertexnum) {
			//顶点映射
			//external idx --> meshlet limit idx map



			std::map<uint, short>& gloidx2idx = gloidx2idxvec[cntid];
			short mapidx = 0;
			for (auto& id : meshlet.externalwireids)
				for (int i = 0; i < Ewires[id].wire.size() - 1; i++)
					gloidx2idx[Ewires[id].wire[i]] = mapidx++;


			int orderofid = 0;
			//对于一条wire只需要存储 wire.size()-1 个 left right corner
			for (auto& id : meshlet.externalwireids) {
				short LorR = NOTINI;

				for (int i = 0; i < Ewires[id].wire.size() - 1; i++) {
					uint startidx = Ewires[id].wire[i];
					uint endidx = Ewires[id].wire[i + 1];
					short idx = FindLeftOrRightIDX_VERSION0(startidx, endidx, meshlet, mesh, gloidx2idx, LorR);
					if (lrref[orderofid] == NOTINI)lrref[orderofid] = LorR;
					else if (lrref[orderofid] != LorR)std::cout << "error in find left or Right" << std::endl;

					if (LorR == RIGHT)Ewires[id].right.push_back(idx);
					else Ewires[id].left.push_back(idx);
				}
				if (LorR == RIGHT)Ewires[id].right.push_back(EXTERNALWIREEND);
				else Ewires[id].left.push_back(EXTERNALWIREEND);
				orderofid++;
			}

		}
		else if ((boundvset.size() + 1) == meshlet.vertexnum) {
			//从startidx开始build loop
			//看论文吧


			for (auto& ver : meshlet.vertexs)
				if (boundvset.find(ver) == boundvset.end()) {
					meshlet.interwire.wire.push_back(ver);
					break;
				}
			//map映射
			std::map<uint, short>& gloidx2idx = gloidx2idxvec[cntid];
			short mapidx = 0;
			gloidx2idx[meshlet.interwire.wire[0]] = mapidx++;
			for (auto& id : meshlet.externalwireids)
				for (int i = 0; i < Ewires[id].wire.size() - 1; i++)
					gloidx2idx[Ewires[id].wire[i]] = mapidx++;


			int orderofid = 0;
			for (auto& id : meshlet.externalwireids) {
				short LorR = NOTINI;

				for (int i = 0; i < Ewires[id].wire.size() - 1; i++) {
					uint startidx = Ewires[id].wire[i];
					uint endidx = Ewires[id].wire[i + 1];
					short idx = FindLeftOrRightIDX_VERSION0(startidx, endidx, meshlet, mesh, gloidx2idx, LorR);
					if (lrref[orderofid] == NOTINI)lrref[orderofid] = LorR;
					else if (lrref[orderofid] != LorR)std::cout << "error in find left or Right" << std::endl;

					if (LorR == RIGHT)Ewires[id].right.push_back(idx);
					else Ewires[id].left.push_back(idx);
				}
				if (LorR == RIGHT)Ewires[id].right.push_back(EXTERNALWIREEND);
				else Ewires[id].left.push_back(EXTERNALWIREEND);
				orderofid++;
			}
		}
		else {
			//find start cycle vertex

			//first elem vertexidx second elem halfedge idx

			//根据论文生成的并不是一个环，而是一条线
			//interlacewire(interwire, mesh, meshlet, boundvset);

			LevelWireGenerator(meshlet, mesh, boundvset);
			std::cout << ' ' << cntid << std::endl;


			std::map<uint, short>& gloidx2idx = gloidx2idxvec[cntid];
			short mapidx = 0;
			for (auto& elem : meshlet.interwire.wire)
				if (elem != NEXTWIRE)
					gloidx2idx[elem] = mapidx++;

			for (auto& id : meshlet.externalwireids)
				for (int i = 0; i < Ewires[id].wire.size(); i++)
					if(gloidx2idx.find(Ewires[id].wire[i])==gloidx2idx.end())
						gloidx2idx[Ewires[id].wire[i]] = mapidx++;

			if (gloidx2idx.size() != meshlet.vertexnum) {
				std::cout << "eroor in map" << std::endl;
				for (auto& elem : meshlet.vertexs)
					if (gloidx2idx.find(elem) == gloidx2idx.end())
						std::cout << elem;
			}
			std::set<uint> facestoremove = meshlet.faces;



			//先build external wire的
			int orderofid = 0;
			for (auto& id : meshlet.externalwireids) {
				short LorR = NOTINI;
				if ((cntid == 1)&&(orderofid==1))
					std::cout << "debug" << std::endl;
				for (int i = 0; i < Ewires[id].wire.size() - 1; i++) {

					uint startidx = Ewires[id].wire[i];
					uint endidx = Ewires[id].wire[i + 1];
	
					short idx = FindLeftOrRightIDX_VERSION1(startidx, endidx, meshlet, mesh, gloidx2idx, LorR,facestoremove);
					if (lrref[orderofid] == NOTINI)lrref[orderofid] = LorR;
					else if (lrref[orderofid] != LorR)std::cout << "error in find left or Right" << std::endl;
					if (LorR == RIGHT) {
						Ewires[id].right.push_back(idx);
					}
					else {
						Ewires[id].left.push_back(idx);
					}
				}
				if (LorR == RIGHT)Ewires[id].right.push_back(EXTERNALWIREEND);
				else Ewires[id].left.push_back(EXTERNALWIREEND);
				orderofid++;
			}

			std::vector<uint> wireref = meshlet.interwire.wire;
			meshlet.interwire.left.resize(wireref.size());
			meshlet.interwire.right.resize(wireref.size());
			for (int wireid = 0; wireid < wireref.size(); wireid++) {
				if (wireref[wireid] == NEXTWIRE)continue;
				if (wireref[wireid + 1] == NEXTWIRE)continue;

				// wireref [wireid] wireref[wireid+1] 为一条线 需要处理
				short rightidx, leftidx;
				FindLeftAndRightIDX_VERSION2(wireref[wireid], wireref[wireid + 1], meshlet, mesh, gloidx2idx, leftidx, rightidx, facestoremove);
				meshlet.interwire.left[wireid] = leftidx;
				meshlet.interwire.right[wireid] = rightidx;
			}

			if (!facestoremove.empty()) {
				meshlet.irregular.push_back(facestoremove.size());
				for (auto& faceid : facestoremove) {
					auto fh = mesh.face_handle(faceid);
					for (const auto& vh : mesh.fv_ccw_range(fh))
						meshlet.irregular.push_back(gloidx2idx[vh.idx()]);
				}
			}
		}
	}

	PackIntoGPUSimple(mesh);

}

LaceWireGenerator::LaceWireGenerator()
{
}


//这个version只适用于没有internal pnt的情况
short LaceWireGenerator::FindLeftOrRightIDX_VERSION0(uint start, uint end, const LaceWire_meshlet& meshlet, const MyMesh& mesh, std::map<uint, short>& idxmap, short& leftorRight)
{
	OpenMesh::VertexHandle vh0(start);


	for (auto it = mesh.cvoh_begin(vh0); it != mesh.cvoh_end(vh0); ++it)
		//找到一条边界out half edge
		if (it->to().idx() == end) {
			int faceid = it->face().idx();
			//假设在右边
			if (meshlet.faces.find(faceid) != meshlet.faces.end()) {
				uint idx = it->next().to().idx();
				leftorRight = LEFT;
				return idxmap[idx];
			}
			else {
			//face在wire左边
				faceid = it->opp().face().idx();
				if (meshlet.faces.find(faceid) == meshlet.faces.end())std::cout << "error" << std::endl;
				uint idx = it->opp().next().to().idx();
				leftorRight = RIGHT;
				return idxmap[idx];
			}
			break;
		}
	std::cout << "error in find left or right version 0" << std::endl;
}

short LaceWireGenerator::FindLeftOrRightIDX_VERSION1(uint start, uint end, const LaceWire_meshlet& meshlet, const MyMesh& mesh, std::map<uint, short>& idxmap, short& leftorRight, std::set<uint>& faces)
{
	OpenMesh::VertexHandle vh0(start);


	for (auto it = mesh.cvoh_begin(vh0); it != mesh.cvoh_end(vh0); ++it)
		//找到一条边界out half edge
		if (it->to().idx() == end) {
			int faceid = it->face().idx();
			//假设在右边
			if (meshlet.faces.find(faceid) != meshlet.faces.end()) {
				uint idx = it->next().to().idx();
				leftorRight = LEFT;
				faces.erase(uint(faceid));
				return idxmap[idx];
			}
			else {
				//face在wire左边
				faceid = it->opp().face().idx();
				if (meshlet.faces.find(faceid) == meshlet.faces.end())std::cout << "error" << std::endl;
				uint idx = it->opp().next().to().idx();
				leftorRight = RIGHT;
				faces.erase(uint(faceid));

				return idxmap[idx];
			}
			break;
		}
	std::cout << "error in find left or right version 0" << std::endl;
}

void LaceWireGenerator::FindLeftAndRightIDX_VERSION2(uint start, uint end, const LaceWire_meshlet& meshlet, const MyMesh& mesh, std::map<uint, short>& idxmap, short& leftidx, short& Rightidx, std::set<uint>& faces)
{
	OpenMesh::VertexHandle vh0(start);
	for (auto it = mesh.cvoh_begin(vh0); it != mesh.cvoh_end(vh0); ++it)
		//找到一条边界out half edge
		if (it->to().idx() == end) {
			int faceid = it->face().idx();
			//假设在右边
			if (meshlet.faces.find(faceid) == meshlet.faces.end()) 
				std::cout << "error in VERSION2" << std::endl;
			uint idx = it->next().to().idx();
			faces.erase(uint(faceid));
			leftidx = idxmap[idx];

			//face在wire左边
			faceid = it->opp().face().idx();
			if (meshlet.faces.find(faceid) == meshlet.faces.end())
				std::cout << "error in VERSION2" << std::endl;
			idx = it->opp().next().to().idx();
			faces.erase(uint(faceid));
			Rightidx = idxmap[idx];
			break;
		}
	return;
}

void LaceWireGenerator::LevelWireGenerator( LaceWire_meshlet& meshlet, const MyMesh& mesh, const std::unordered_set<uint>& boundset)
{
	std::vector<std::set<uint>> levelofPnts;

	LevelOfontGeneratr(levelofPnts,meshlet, mesh, boundset);
	std::vector<std::deque<uint>> internalwires;

	int remainnums = meshlet.vertexnum - boundset.size();
	while (remainnums != 0) {
		std::deque<uint> onewire;
		uint startidx = GetstartPnts(levelofPnts);
		remainnums--;
		onewire.push_front(startidx);
		GenerateOneWire(onewire, remainnums, levelofPnts, mesh);
		internalwires.push_back(onewire);
	}

	for (auto& wire : internalwires) {
		for (auto& elem : wire)
			meshlet.interwire.wire.push_back(elem);
		meshlet.interwire.wire.push_back(NEXTWIRE);
	}

	int cnt = boundset.size();
	for (auto& wire : internalwires)
		cnt += wire.size();
	if (cnt != meshlet.vertexnum)
		std::cout << "Error in internal wire generta" << std::endl;


	std::cout << "line num of meshlet id   " << internalwires.size();
}

void LaceWireGenerator::LevelOfontGeneratr(std::vector<std::set<uint>>& levelpnts, const LaceWire_meshlet& meshlet, const MyMesh& mesh, const std::unordered_set<uint>& boundset)
{
	std::set<uint> elems;
	std::set<uint> remainpnts{meshlet.vertexs};
	for (auto& elem : boundset) {
		remainpnts.erase(elem);
		elems.insert(elem);
	}

	while (!remainpnts.empty())
	{
		std::set<uint> pnts;

		for (auto& elem : elems) {
			auto vh = mesh.vertex_handle(elem);
			for (const auto& vvit : mesh.vv_range(vh))
				if (remainpnts.find(vvit.idx()) != remainpnts.end()) {
					pnts.insert(vvit.idx());
					remainpnts.erase(vvit.idx());
				}
		}
		levelpnts.push_back(pnts);
		elems.clear();
		elems = pnts;
	}
}

uint LaceWireGenerator::GetstartPnts(std::vector<std::set<uint>>& levelpnts)
{
	for (auto& elem : levelpnts)
		if (!elem.empty()) {
			auto pnt = *elem.begin();
			elem.erase(pnt);
			return pnt;
		}
}

void LaceWireGenerator::GenerateOneWire(std::deque<uint>& onewire, int& remainnums, std::vector<std::set<uint>>& levelofPnts, const MyMesh& mesh)
{
	while (remainnums!=0)
	{
		bool flag = false;
		uint startidx = onewire.front();
		auto vh = mesh.vertex_handle(startidx);
		for (const auto& vvit : mesh.vv_cw_range(vh))
			if (IsNextPnt(vvit.idx(), levelofPnts)) {
				remainnums--;
				onewire.push_front(vvit.idx());
				flag = true;
				break;
			}

		if (flag)
			continue;

		uint endidx = onewire.back();
		vh = mesh.vertex_handle(endidx);
		for(const auto&vvit:mesh.vv_range(vh))
			if (IsNextPnt(vvit.idx(), levelofPnts)) {
				remainnums--;
				onewire.push_back(vvit.idx());
				flag = true;
				break;
			}
		if (flag)
			continue;

		//两端都插入失败
		return;
	}
	return;
}

bool LaceWireGenerator::IsNextPnt(uint vidx, std::vector<std::set<uint>>& levelofPnts)
{
	for(auto&level:levelofPnts)
		if ((!level.empty())&&(level.find(vidx) != level.end())) {
			level.erase(vidx);
			return true;
		}
	return false;
}

void LaceWireGenerator::PackIntoGPUSimple(const MyMesh& mesh)
{
	uint vertexbegin = 0;
	uint primbegin = 0;
	int cntid = -1;
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<float> dis(0.0, 1.0);

	for (const auto& meshlet : meshlets) {
		//geoinfo
		cntid++;
		uint vertexcnt = 0;

		std::map<uint, short>& mapref = gloidx2idxvec[cntid];
		for (const auto& elem : meshlet.interwire.wire) {
			if (elem == NEXTWIRE)continue;
			const auto& vh = mesh.vertex_handle(elem);
			auto pnt = mesh.point(vh);
			geoinfo.emplace_back(pnt[0], pnt[1], pnt[2],1.0f);
			vertexcnt++;
		}

		for (const auto& id : meshlet.externalwireids) {
			for (const auto& elem : Ewires[id].wire) {
				const auto& vh = mesh.vertex_handle(elem);
				auto pnt = mesh.point(vh);
				geoinfo.emplace_back(pnt[0], pnt[1], pnt[2], 1.0f);
				vertexcnt++;
			}
			geoinfo.pop_back();
			vertexcnt--;
		}

		if (vertexcnt != meshlet.vertexnum)
			std::cout << "error in pack geoinfo" << std::endl;
		
		for (uint i = 0; i < meshlet.interwire.left.size() - 1; i++) {
			if (meshlet.interwire.wire[i] == NEXTWIRE)continue;
			if (meshlet.interwire.wire[i + 1] == NEXTWIRE)continue;
			uint leftidx = meshlet.interwire.left[i];
			uint rightidx = meshlet.interwire.right[i];
			uint elem = meshlet.interwire.wire[i];
			uint elemadd = meshlet.interwire.wire[i + 1];
			priminfo.push_back((unsigned char)mapref[elem]);
			priminfo.push_back((unsigned char)mapref[elemadd]);
			priminfo.push_back((unsigned char)leftidx);

			priminfo.push_back((unsigned char)mapref[elemadd]);
			priminfo.push_back((unsigned char)mapref[elem]);
			priminfo.push_back((unsigned char)rightidx);
		}


		short irsize = meshlet.irregular[0];
		for (int i = 0; i < irsize; i++) {
			priminfo.push_back((unsigned char)meshlet.irregular[3 * i + 1]);
			priminfo.push_back((unsigned char)meshlet.irregular[3 * i + 2]);
			priminfo.push_back((unsigned char)meshlet.irregular[3 * i + 3]);

		}

		for (int i = 0; i < meshlet.externalwireids.size(); i++) {
			short sign = meshlet.externalLR[i];
			int eid = meshlet.externalwireids[i];
			for (int j = 0; j < Ewires[eid].wire.size() - 1; j++) {
				if (sign == LEFT) {
					unsigned char leftidx = Ewires[eid].left[j];
					uint elem = Ewires[eid].wire[j];
					uint elemadd = Ewires[eid].wire[j+1];
					priminfo.push_back((unsigned char)mapref[elem]);
					priminfo.push_back((unsigned char)mapref[elemadd]);
					priminfo.push_back(leftidx);
				}
				else {
					unsigned char rightidx = Ewires[eid].right[j];
					uint elem = Ewires[eid].wire[j];
					uint elemadd = Ewires[eid].wire[j+1];
					priminfo.push_back((unsigned char)mapref[elemadd]);
					priminfo.push_back((unsigned char)mapref[elem]);
					priminfo.push_back(rightidx);
				}
			}
		}



		Simple_meshlet temp;
		temp.vertex_begin = vertexbegin;
		temp.vertex_cnt = vertexcnt;


		temp.primbegin = primbegin;
		temp.primcnt = (priminfo.size() - primbegin) / 3;
		temp.color = vec4(dis(gen), dis(gen), dis(gen), 1.0f);

		vertexbegin = vertexbegin + vertexcnt;
		primbegin = priminfo.size();
		simplemeshlets.push_back(temp);
	}

}

//void LaceWireGenerator::Interlacewire(std::vector<std::array<uint, 2>>& wirerecord, const MyMesh& mesh, LaceWire_meshlet& meshlet, const std::unordered_set<uint>& boundset)
//{
//	uint starthlidx;
//	uint startidx = FindStartVertex(meshlet, mesh, boundset, starthlidx);
//	uint nextidx, nexthlidx;
//	uint intervertexnum = meshlet.vertexnum - boundset.size();
//	wirerecord.push_back(std::array<uint, 2>{ startidx, 0});
//	meshlet.interwire.wire.push_back(startidx);
//	std::set<uint> alreadyin{ startidx };
//
//	while (Internextpntsearch(mesh,boundset,startidx,starthlidx,nextidx,nexthlidx,alreadyin)) {
//		if (meshlet.vertexs.find(nextidx) == meshlet.vertexs.end())
//			std::cout << "error in occur outside pnts";
//
//		auto edgehl = mesh.halfedge_handle(nexthlidx);
//		edgehl = mesh.prev_halfedge_handle(edgehl);
//		edgehl = mesh.opposite_halfedge_handle(edgehl);
//		if (mesh.to_vertex_handle(edgehl).idx() != nextidx)
//			std::cout << "error" << std::endl;
//		if(mesh.from_vertex_handle(edgehl).idx()!=startidx)
//			std::cout << "error" << std::endl;
//		std::array<uint, 2> temp;
//		temp[0] = nextidx; temp[1] = edgehl.idx();
//
//
//		wirerecord.push_back(temp);
//		alreadyin.insert(nextidx);
//		meshlet.interwire.wire.push_back(nextidx);
//
//		startidx = nextidx;
//		starthlidx = nexthlidx;
//
//	}
//
//	if ((alreadyin.size() + boundset.size()) != meshlet.vertexnum) {
//		std::cout << "not a full ring" << std::endl;
//	}
//
//}
//
//bool LaceWireGenerator::Internextpntsearch(const MyMesh& mesh, const std::unordered_set<uint>& boundset, uint& startidx, uint& starthlidx, uint& nextidx, uint& nexthlidx, const std::set<uint>& pre)
//{
//	auto hl = mesh.halfedge_handle(starthlidx);
//	
//	do {
//		hl = mesh.next_halfedge_handle(hl);
//		hl = mesh.opposite_halfedge_handle(hl);
//		hl = mesh.next_halfedge_handle(hl);
//		auto vh = mesh.to_vertex_handle(hl);
//		if (	(boundset.find(vh.idx()) == boundset.end())
//				&&//vh.idx既不在boundset又不在preset中间
//			(pre.find(vh.idx())==pre.end())	) {
//			nextidx = vh.idx();
//			nexthlidx = mesh.prev_halfedge_handle(hl).idx();
//			return true;
//		}
//
//
//	} while (hl.idx() != starthlidx);
//		return false;
//}

//uint LaceWireGenerator::FindStartVertex(const LaceWire_meshlet& meshlet, const MyMesh& mesh, const std::unordered_set<uint>& boundset, uint& hlidx) {
//	//第二次时出现错误
//	for (auto& id : meshlet.externalwireids)
//		for (int i = 0; i < Ewires[i].wire.size() - 1; i++) {
//			int startvid = Ewires[i].wire[i];
//			int endvid = Ewires[i].wire[i + 1];
//
//			OpenMesh::VertexHandle vh0(startvid);
//
//			OpenMesh::VertexHandle vh1(endvid);
//
//			for (auto it = mesh.cvoh_begin(vh0); it != mesh.cvoh_end(vh0); ++it)
//				//找到一条边界out half edge
//				if (it->to().idx() == endvid) {
//					int faceid = it->face().idx();
//
//					//假设在start -->  end的右边
//					if (meshlet.faces.find(faceid) != meshlet.faces.end()) {
//						uint vidx = it->next().to().idx();
//						if (boundset.find(vidx) == boundset.end()) {
//							hlidx = it->idx();
//							return vidx;
//						}
//						else
//							continue;
//					}
//					//在左边
//					faceid = it->opp().face().idx();
//					if (meshlet.faces.find(faceid) == meshlet.faces.end())
//						std::cout << "error" << std::endl;
//					uint vidx = it->opp().next().to().idx();
//					if (boundset.find(vidx) == boundset.end()) {
//						hlidx = it->opp().idx();
//						return vidx;
//					}
//					else
//						continue;
//				}
//
//		}
//	return NOINTERNALVERTEX;
//
//}