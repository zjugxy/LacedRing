#include "NewLWGenerator.h"
#define EMPTYWIRE 0xFFFFFFFF
#define WIREEND 0xFF
#define REVERSEFLAG 0x80000000

#include<random>

NewLWGenerator::NewLWGenerator(const NewCluster& nclu)
{
	nummeshlet = nclu.mymeshlets.size();
	targets.resize(nummeshlet);
	internalbuilders.resize(nummeshlet);

	for (uint i = 0; i < nummeshlet; ++i) {
		internalbuilders[i].borderset = nclu.buildsets[i].borderset;
		internalbuilders[i].cornerset = nclu.buildsets[i].cornerset;
		internalbuilders[i].faces = nclu.mymeshlets[i].faces;
		internalbuilders[i].vertices = nclu.mymeshlets[i].vertices;

		targets[i].ewireidx = nclu.buildsets[i].ewiresid;
		targets[i].reverse = nclu.buildsets[i].ewirereverse;
	}

	gloEwires.resize(nclu.gloewires.size());
	for (uint i = 0; i < nclu.gloewires.size(); ++i) {
		gloEwires[i].vertex = nclu.gloewires[i];
		gloEwires[i].left.resize(gloEwires[i].vertex.size());
		gloEwires[i].right.resize(gloEwires[i].vertex.size());
	}

}

void NewLWGenerator::FUNC(const MyMesh& mesh)
{
	InternalWireGenerator(mesh);
	LeftRightGenertor(mesh);
	PackSimpleLaceWire(mesh);
}

void NewLWGenerator::InternalWireGenerator(const MyMesh& mesh)
{
	for (uint mid = 0; mid < nummeshlet; ++mid) {
		auto& builder = internalbuilders[mid];
		auto& target = targets[mid];
		//没有内在顶点
		if (builder.borderset.size() == builder.vertices.size()) {
			target.InternalLW.vertex.push_back(EMPTYWIRE);
			target.InternalLW.left.push_back(WIREEND);
			target.InternalLW.right.push_back(WIREEND);
			continue;
		}
		//只有一个内在顶点
		else if (builder.borderset.size() == (builder.vertices.size() - 1)) {
			for(auto&ver:builder.vertices)
				if (builder.borderset.find(ver) == builder.borderset.end()) {
					target.InternalLW.vertex.push_back(ver);
					break;
				}
			assert(target.InternalLW.vertex.size() == 1);
			target.InternalLW.left.push_back(WIREEND);
			target.InternalLW.right.push_back(WIREEND);
			continue;
		}
		else {
			//建立interwire
			std::vector<std::unordered_set<uint >> levelofpnts;
			std::unordered_set<uint> pntstoremove;
			LevelPntGen(builder.borderset, builder.vertices, levelofpnts, mesh);

			for (auto& level : levelofpnts)
				for (auto& ver : level)
					pntstoremove.insert(ver);

			while (!pntstoremove.empty())
			{
				std::deque<uint> singlewire;
				for (auto& level : levelofpnts)
					if (!level.empty()) {
						singlewire.push_back(*level.begin());
						level.erase(singlewire.back());
						pntstoremove.erase(singlewire.back());
						break;
					}

				PushAtBack(levelofpnts, mesh, singlewire, pntstoremove);
				PushAtFront(levelofpnts, mesh, singlewire, pntstoremove);
				for (auto& ver : singlewire)
					target.InternalLW.vertex.push_back(ver);
				target.InternalLW.left.resize(target.InternalLW.vertex.size());
				target.InternalLW.right.resize(target.InternalLW.vertex.size());
				target.InternalLW.left.back() = WIREEND;
				target.InternalLW.right.back() = WIREEND;
			}
		}
	}
}

void NewLWGenerator::LeftRightGenertor(const MyMesh& mesh)
{
	vidxmaps.resize(nummeshlet);

	for (uint mid = 0; mid < nummeshlet; ++mid) {
		std::unordered_map<uint, uchar>& vidxmap = vidxmaps[mid];
		auto& builder = internalbuilders[mid];
		auto& target = targets[mid];
		auto& intermesh = target.InternalLW;

		auto facetoremove = builder.faces;
		//填充vidxmap
		{
			uchar vidx = 0;
			for (auto& ver : target.InternalLW.vertex)
				if (ver == EMPTYWIRE)
					break;
				else
					vidxmap[ver] = vidx++;

			for (uint i = 0; i < target.ewireidx.size(); ++i) {
				auto& evertex = gloEwires[target.ewireidx[i]].vertex;
				if (target.reverse[i] == false) {
					for (uint j = 0; j < evertex.size() - 1; ++j) {
						if (vidxmap.find(evertex[j]) == vidxmap.end())
							vidxmap[evertex[j]] = vidx++;
						else
							std::cout << "error may occur in vidxmap sharp"<<mid;
					}
				}
				else {
					for (uint j = evertex.size() - 1; j != 0; --j) {
						if (vidxmap.find(evertex[j]) == vidxmap.end())
							vidxmap[evertex[j]] = vidx++;
						else
							std::cout << "error may occur in vidxmap sharp"<<mid;
					}
				}
			}
			assert(vidxmap.size() == internalbuilders[mid].vertices.size());
		}

		//开始填充internal wire的left right
		{
			if (intermesh.vertex[0] == EMPTYWIRE) {
				assert(intermesh.vertex.size() == 1);
				intermesh.left[0] = WIREEND;
				intermesh.right[0] = WIREEND;
			}
			else {
				for (uint i = 0; i < intermesh.vertex.size(); ++i) {
					if (intermesh.left[i] == WIREEND)continue;

					uint startidx = intermesh.vertex[i];
					uint endidx = intermesh.vertex[i+1];

					auto vh = mesh.vertex_handle(startidx);
					for (auto vohit = mesh.cvoh_iter(vh); vohit.is_valid(); ++vohit) {
						auto toidx = vohit->to().idx();
						if (toidx != endidx)
							continue;
						else {
							intermesh.left[i] = vidxmap[vohit->next().to().idx()];
							intermesh.right[i] = vidxmap[vohit->opp().next().to().idx()];
							facetoremove.erase(vohit->face().idx());
							facetoremove.erase(vohit->opp().face().idx());
							break;
						}
					}
				}
			}
		}

		//开始填充external wire

		for (uint i = 0; i < target.ewireidx.size(); ++i) {
			auto& evertex = gloEwires[target.ewireidx[i]].vertex;
			auto& emesh = gloEwires[target.ewireidx[i]];
			if (target.reverse[i] == false) {
				//顺序遍历

				for (uint j = 0; j < evertex.size() - 1; ++j) {
					uint startidx = emesh.vertex[j];
					uint endidx = emesh.vertex[j + 1];

					auto vh = mesh.vertex_handle(startidx);
					for (auto vohit = mesh.cvoh_iter(vh); vohit.is_valid(); ++vohit) {
						auto toidx = vohit->to().idx();
						if (toidx != endidx)
							continue;
						else {
							auto vidx = vohit->next().to().idx();
							assert(vidxmap.find(vidx) != vidxmap.end());
							emesh.left[j] = vidxmap[vohit->next().to().idx()];
							facetoremove.erase(vohit->face().idx());
							break;
						}
					}
				}
				emesh.left[evertex.size() - 1] = WIREEND;
			}
			else {
				//逆序遍历
				for (uint j = evertex.size()-1; j != 0; --j) {
					uint startidx = emesh.vertex[j];
					uint endidx = emesh.vertex[j - 1];

					auto vh = mesh.vertex_handle(startidx);
					for (auto vohit = mesh.cvoh_iter(vh); vohit.is_valid(); ++vohit) {
						auto toidx = vohit->to().idx();
						if (toidx != endidx)
							continue;
						else {
							auto vidx = vohit->next().to().idx();
							assert(vidxmap.find(vidx) != vidxmap.end());
							emesh.right[j] = vidxmap[vohit->next().to().idx()];
							facetoremove.erase(vohit->face().idx());

							break;
						}
					}
				}
				emesh.right[0] = WIREEND;
			}
		}

		for (auto& face : facetoremove) {
			auto fh = mesh.face_handle(face);
			for (auto fvit = mesh.cfv_iter(fh); fvit.is_valid(); ++fvit) {
				intermesh.irregular.push_back(vidxmap[fvit->idx()]);
			}
		}

	}
}

void NewLWGenerator::LevelPntGen(const std::unordered_set<uint>& bordset, const std::unordered_set<uint>& verset,
	std::vector<std::unordered_set<uint>>& lpnts, const MyMesh& mesh){
	std::unordered_set<uint> remainpnts(verset);
	std::unordered_set<uint> elems;

	for (const auto& elem : bordset) {
		remainpnts.erase(elem);
		elems.insert(elem);
	}

	while (!remainpnts.empty())
	{
		std::unordered_set<uint> pnts;

		for (auto& elem : elems) {
			auto vh = mesh.vertex_handle(elem);
			for (const auto& vvit : mesh.vv_range(vh))
				if (remainpnts.find(vvit.idx()) != remainpnts.end()) {
					pnts.insert(vvit.idx());
					remainpnts.erase(vvit.idx());
				}
		}
		lpnts.push_back(pnts);
		elems.clear();
		elems = pnts;
	}
}

void NewLWGenerator::PushAtBack(std::vector<std::unordered_set<uint>>& lpnts, const MyMesh& mesh, std::deque<uint>& singlewire, std::unordered_set<uint>& pntstoremove)
{
	do {
		auto vidx = singlewire.back();
		auto vh = mesh.vertex_handle(vidx);

		//找顺势正遍历的最后一个wo'r'k的顶点
		std::vector<int> record;

		auto vvit = mesh.cvv_cwiter(vh);
		for (; vvit.is_valid(); ++vvit) {
			record.push_back(vvit->idx());
		}

		int nextpnt = -1;
		int startidx = -1;
		//需要在removepnts中间开始遍历，
		for (int i = 0; i < record.size(); i++)
			if (pntstoremove.find(record[i]) == pntstoremove.end())
				startidx = i;
		assert(startidx != -1);
		//从 startidx开始
		for (int i = 1; i < record.size(); i++) {
			auto ver = record[(i + startidx) % record.size()];
			if (pntstoremove.find(ver) != pntstoremove.end())
				nextpnt = ver;
		}

		if (nextpnt == -1)
			break;
		else {
			singlewire.push_back(nextpnt);
			pntstoremove.erase(nextpnt);
			for (auto& level : lpnts)
				level.erase(nextpnt);
		}
	} while (!pntstoremove.empty());

}

void NewLWGenerator::PushAtFront(std::vector<std::unordered_set<uint>>& lpnts, const MyMesh& mesh, std::deque<uint>& singlewire, std::unordered_set<uint>& pntstoremove)
{

	do {
		auto vidx = singlewire.front();
		auto vh = mesh.vertex_handle(vidx);

		//找逆时针正遍历的最后一个wo'r'k的顶点
		std::vector<int> record;

		auto vvit = mesh.cvv_ccwiter(vh);
		for (; vvit.is_valid(); ++vvit) {
			record.push_back(vvit->idx());
		}

		int nextpnt = -1;
		int startidx = -1;
		//需要在removepnts中间开始遍历，
		for (int i = 0; i < record.size(); i++)
			if (pntstoremove.find(record[i]) == pntstoremove.end())
				startidx = i;
		assert(startidx != -1);
		//从 startidx开始
		for (int i = 1; i < record.size(); i++) {
			auto ver = record[(i + startidx) % record.size()];
			if (pntstoremove.find(ver) != pntstoremove.end())
				nextpnt = ver;
		}

		if (nextpnt == -1)
			break;
		else {
			singlewire.push_front(nextpnt);
			pntstoremove.erase(nextpnt);
			for (auto& level : lpnts)
				level.erase(nextpnt);
		}
	} while (!pntstoremove.empty());
}

void NewLWGenerator::PackSimpleLaceWire(const MyMesh& mesh)
{

	//
	//geoinfo要使用原始方法pack,那么可能要处理冲突问题（但是只在meshlet为50的meshlet出现了哪个情况
	//环怎么处理，2 loop-->有问题
	//

	simplemeshlets.resize(nummeshlet);
	uint vertexbegin = 0;
	uint primbegin = 0;
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<float> dis(0.0, 1.0);


	for (uint mid = 0; mid < nummeshlet; ++mid) {

		auto& vidxmap = vidxmaps[mid];
		auto& meshlet = targets[mid];
		auto& intermeshlet = meshlet.InternalLW;

		geoinfo.resize(vertexbegin + vidxmap.size());
		//for (auto& elem : vidxmap) {
		//	auto vh = mesh.vertex_handle(elem.first);
		//	auto pnt = mesh.point(vh);
		//	geoinfo[vertexbegin + elem.second] = vec4{ pnt[0],pnt[1],pnt[2],1.0 };
		//}

		uint vcnt = 0;
		for (auto& ver : meshlet.InternalLW.vertex) {
			if (ver == EMPTYWIRE)
				break;
			auto vh = mesh.vertex_handle(ver);
			auto pnt = mesh.point(vh);
			geoinfo[vertexbegin + vcnt] = vec4{ pnt[0],pnt[1],pnt[2],1.0 };
			vcnt++;
		}

		for (uint eid = 0; eid < meshlet.ewireidx.size(); ++eid) {
			auto& ewire = gloEwires[meshlet.ewireidx[eid]];
			if (meshlet.reverse[eid] == false) {
				for (uint i = 0; i < ewire.vertex.size() - 1; ++i) {
					auto vh = mesh.vertex_handle(ewire.vertex[i]);
					auto pnt = mesh.point(vh);
					geoinfo[vertexbegin + vcnt] = vec4{ pnt[0],pnt[1],pnt[2],1.0 };
					vcnt++;
				}
			}
			else {
				for (uint i = ewire.vertex.size() - 1; i !=0 ; --i) {
					auto vh = mesh.vertex_handle(ewire.vertex[i]);
					auto pnt = mesh.point(vh);
					geoinfo[vertexbegin + vcnt] = vec4{ pnt[0],pnt[1],pnt[2],1.0 };
					vcnt++;
				}
			}
		}
		assert(vcnt == vidxmap.size());





		for (uint i = 0; i < intermeshlet.vertex.size(); ++i) {
			if (intermeshlet.vertex[i] == EMPTYWIRE)
				break;
			if (intermeshlet.left[i] == WIREEND)
				continue;
			priminfo.push_back(i);
			priminfo.push_back(i + 1);
			priminfo.push_back(intermeshlet.left[i]);

			priminfo.push_back(i + 1);
			priminfo.push_back(i);
			priminfo.push_back(intermeshlet.right[i]);
		}


		uchar idxstart = 0;
		uchar ivertexnum;
		for (auto& elem : intermeshlet.irregular)
			priminfo.push_back(elem);


		if (intermeshlet.vertex[0] == EMPTYWIRE)
			ivertexnum = 0;
		else
			ivertexnum = intermeshlet.vertex.size();

		uchar evertexnum = vidxmap.size() - ivertexnum;


		//if (ivertexnum == 0)
		//	std::cout << "debug" << mid << std::endl;
		//if (mid == 56)
		//	std::cout << "try" << std::endl;


		for (uint eid = 0; eid < meshlet.ewireidx.size(); ++eid) {
			auto& ewire = gloEwires[meshlet.ewireidx[eid]];
			if (meshlet.reverse[eid] == false) {
				for (uint i = 0; i < ewire.vertex.size()-1; ++i) {
					
					priminfo.push_back(ivertexnum + (idxstart + i) % evertexnum);
					priminfo.push_back(ivertexnum + (idxstart + i + 1) % evertexnum);
					priminfo.push_back(ewire.left[i]);

					if ((ivertexnum + (idxstart + i) % evertexnum) != vidxmap[ewire.vertex[i]])
						std::cout << "mid " << mid << " " << (uint)(ivertexnum + (idxstart + i) % evertexnum)
						<< " " << (uint)(vidxmap[ewire.vertex[i]]) << std::endl;
					if ((ivertexnum + (idxstart + i+1) % evertexnum) != vidxmap[ewire.vertex[i+1]])
						std::cout << "mid " << mid << " " << (ivertexnum + (idxstart + i+1) % evertexnum)
						<< " " << (uint)(vidxmap[ewire.vertex[i+1]]) << std::endl;
				}
			}
			else {
				for (uint i = 0; i < ewire.vertex.size() - 1; ++i) {
					priminfo.push_back(ivertexnum + (idxstart + i) % evertexnum);
					priminfo.push_back(ivertexnum + (idxstart + i +1) % evertexnum);
					priminfo.push_back(ewire.right[ewire.vertex.size() -1 -i]);

					if ((ivertexnum + (idxstart + i) % evertexnum) != vidxmap[ewire.vertex[ewire.vertex.size() - 1 - i]])
						std::cout <<"mid "<<mid<< " " << ivertexnum + (idxstart + i) % evertexnum
						<< " "<< uint(vidxmap[ewire.vertex[ewire.vertex.size() - 1 - i]]) << std::endl;
					if ((ivertexnum + (idxstart + i + 1) % evertexnum) != vidxmap[ewire.vertex[ewire.vertex.size() - 1 - i - 1]])
						std::cout << "mid " << mid << " " << ivertexnum + (idxstart + i+1) % evertexnum
						<< " " << uint(vidxmap[ewire.vertex[ewire.vertex.size() - 2 - i]]) << std::endl;
				}
			}
			idxstart += ewire.vertex.size() - 1;
		}


		simplemeshlets[mid].color = vec4(dis(gen), dis(gen), dis(gen), 1.0f);
		simplemeshlets[mid].vertex_begin = vertexbegin;
		simplemeshlets[mid].vertex_cnt = vidxmap.size();

		simplemeshlets[mid].primbegin = primbegin;
		simplemeshlets[mid].primcnt = (priminfo.size() - primbegin) / 3;

		primbegin = priminfo.size();
		vertexbegin += vidxmap.size();

	}

}

void NewLWGenerator::PackGPULW(const MyMesh& mesh)
{
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<float> dis(0.0, 1.0);

	//
	std::vector<uint> interrecord;
	std::vector<uint> exterrecord;
	//pack intermesh
	for (uint mid = 0; mid < targets.size(); mid++) {
		interrecord.push_back(Intermesh.size());

		auto& target = targets[mid];

		std::vector<uchar> vec;
		vec.push_back(static_cast<uchar>(vidxmaps[mid].size()));
		vec.push_back(uchar(0));
		vec.push_back(static_cast<uchar>(target.InternalLW.irregular.size() / 3));

		assert(vidxmaps[mid].size() == target.InternalLW.left.size());
		if (target.InternalLW.vertex[0] != EMPTYWIRE) {
			for (auto& elem : target.InternalLW.left)
				vec.push_back(elem);
			for (auto& elem : target.InternalLW.right)
				vec.push_back(elem);
		}
		for (auto& elem : target.InternalLW.irregular)
			vec.push_back(elem);

		while (vec.size()%4 != 0)
		{
			vec.push_back(0);
		}
		vec[1] = vec.size() / 4;

		for (uint i = 0; i < vec.size(); i += 4) {
			Intermesh.push_back(PackChar4Float(vec[i], vec[i + 1], vec[i + 2], vec[i + 3]));
		}
		Intermesh.push_back(dis(gen));
		Intermesh.push_back(dis(gen));
		Intermesh.push_back(dis(gen));

		for (auto& elem : target.InternalLW.vertex) {
			if (elem == EMPTYWIRE)
				break;
			else {
				auto pnt = mesh.point(mesh.vertex_handle(elem));
				Intermesh.push_back(pnt[0]);
				Intermesh.push_back(pnt[1]);
				Intermesh.push_back(pnt[2]);
			}
		}

	}

	//pack external
	for (auto& ewire : gloEwires) {
		exterrecord.push_back(Extermesh.size());
		std::vector<uchar> vec;
		vec.push_back(static_cast<uchar>(ewire.vertex.size()));
		vec.push_back(uchar(0));

		for (auto& elem : ewire.left)
			vec.push_back(elem);
		for (auto& elem : ewire.right)
			vec.push_back(elem);
		while (vec.size() % 4 != 0)
		{
			vec.push_back(0);
		}
		vec[1] = vec.size() / 4;
		for (uint i = 0; i < vec.size(); i += 4) {
			Extermesh.push_back(PackChar4Float(vec[i], vec[i + 1], vec[i + 2], vec[i + 3]));
		}

		for (auto& elem : ewire.vertex) {
			auto pnt = mesh.point(mesh.vertex_handle(elem));
			Intermesh.push_back(pnt[0]);
			Intermesh.push_back(pnt[1]);
			Intermesh.push_back(pnt[2]);
		}
	}
	
	assert(Extermesh.size() < 0xEFFFFFFF);

	for (uint mid = 0; mid < targets.size(); ++mid) {
		PositionInfo.push_back(DesInfo.size());
		uchar vertexcnt = static_cast<uchar>(vidxmaps[mid].size());
		uchar ewirecnt = static_cast<uchar>(targets[mid].ewireidx.size());
		uchar primcnt = 2 * (targets[mid].InternalLW.vertex.size() - 1) + vertexcnt - targets[mid].InternalLW.vertex.size();
		DesInfo.push_back(PackChar4Uint(vertexcnt, ewirecnt, primcnt, 0));

		DesInfo.push_back(interrecord[mid]);

		for (uint i = 0; i < targets[mid].ewireidx.size(); i++) {
			if (targets[mid].reverse[i] == false) {
				DesInfo.push_back(exterrecord[targets[mid].ewireidx[i]]);
			}
			else {
				DesInfo.push_back(exterrecord[targets[mid].ewireidx[i]]|REVERSEFLAG);
			}
		}
	}
}

std::array<uchar, 4> NewLWGenerator::DiscomposeUint(uint value) {
	std::array<uchar, 4> result;
	result[0] = static_cast<uchar>((value >> 24) & 0xFF);
	result[1] = static_cast<uchar>((value >> 16) & 0xFF);
	result[2] = static_cast<uchar>((value >> 8) & 0xFF);
	result[3] = static_cast<uchar>(value & 0xFF);
	return result;
}

float NewLWGenerator::PackChar4Float(uchar c0, uchar c1, uchar c2, uchar c3) {
	uint32_t packedValue = (static_cast<uint32_t>(c0) << 24) |
		(static_cast<uint32_t>(c1) << 16) |
		(static_cast<uint32_t>(c2) << 8) |
		static_cast<uint32_t>(c3);
	return reinterpret_cast<float&>(packedValue);
}

uint NewLWGenerator::PackChar4Uint(uchar c0, uchar c1, uchar c2, uchar c3) {
	uint32_t packedValue = (static_cast<uint32_t>(c0) << 24) |
		(static_cast<uint32_t>(c1) << 16) |
		(static_cast<uint32_t>(c2) << 8) |
		static_cast<uint32_t>(c3);
	return reinterpret_cast<uint&>(packedValue);
}	