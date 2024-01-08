#include "NewLWGenerator.h"
#define EMPTYWIRE 0xFFFFFFFF
#define WIREEND 0xFF
#define REVERSEFLAG 0x80000000
#include<fstream>
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

	//VertexQuantization(mesh);
	NewVertexQuantization(mesh);
	PackGPULW(mesh);
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

void NewLWGenerator::VertexQuantization(const MyMesh& mesh)
{
	//compute full presion
	packexter.resize(gloEwires.size());
	packinter.resize(targets.size());

	bool firstbox = true;
	UnitBox MeshScaleBox;
	UnitBox MeshTransBox;


	for (uint eid = 0; eid < gloEwires.size(); ++eid) {
		const std::vector<uint>& vertexset = gloEwires[eid].vertex;

		if (vertexset.size() < 3) {
			packexter[eid].needtransform = false;
			for (auto ver : vertexset) {
				auto vh = mesh.vertex_handle(ver);
				auto pnt = mesh.point(vh);
				packexter[eid].nopackgeo.emplace_back(pnt[0], pnt[1], pnt[2]);
			}
			continue;
		}
		assert(vertexset.size() != 1);

		Eigen::MatrixXd WireVertices(3, vertexset.size());
		for (uint i = 0; i < vertexset.size(); ++i) {
			auto ver = vertexset[i];
			auto vh = mesh.vertex_handle(ver);
			auto pnt = mesh.point(vh);
			WireVertices.col(i) << pnt[0], pnt[1], pnt[2];
		}

		Eigen::Vector3d centroid = WireVertices.rowwise().mean().transpose();
		Eigen::MatrixXd centeredPoints = WireVertices.colwise() - centroid;
		Eigen::JacobiSVD<Eigen::MatrixXd> svd(centeredPoints, Eigen::ComputeThinU | Eigen::ComputeThinV);
		Eigen::Vector3d normal = svd.matrixU().col(2).normalized();

		Eigen::Vector3d starty, startz;
		getStartAxis(starty, startz, normal);

		Box2d globox{ -1000000.0,-1000000.0,1000000.0,1000000.0 };
		uint glocnt = 0;
		float anglediff = 15.0f;

		for (uint anglecnt = 0; anglecnt < uint(90.0f / anglediff); anglecnt++) {
			float angle = anglecnt * anglediff;
			float radian = angle / 180.0f * PI;
			Eigen::Vector3d yaxis = cos(radian) * starty + sin(radian) * startz;
			Eigen::Vector3d zaxis = cos(radian) * startz - sin(radian) * starty;
			assert(yaxis.dot(zaxis) < 0.0001f);

			//计算2维包围盒
			Box2d box;

			for (auto ver : vertexset) {
				auto vh = mesh.vertex_handle(ver);
				auto pnt = mesh.point(vh);
				Eigen::Vector3d temppnt{ pnt[0],pnt[1],pnt[2] };
				Eigen::Vector3d p2p = temppnt - centroid;
				float yvalue = p2p.dot(yaxis);
				float zvalue = p2p.dot(zaxis);
				box.refresh(yvalue, zvalue);
			}
			if (box.getSumLength() < globox.getSumLength()) {
				globox = box;
				glocnt = anglecnt;
			}
		}

		Eigen::Vector3d finalxaxis = normal.normalized();
		Eigen::Vector3d finalyaxis = cos(glocnt * anglediff / 180.0 * PI) * starty + sin(glocnt * anglediff / 180.0 * PI) * startz;
		Eigen::Vector3d finalzaxis = cos(glocnt * anglediff / 180.0 * PI) * startz - sin(glocnt * anglediff / 180.0 * PI) * starty;
		Eigen::Vector3d newcentroid = centroid + (globox.maxy + globox.miny) / 2.0f * finalyaxis + (globox.maxz + globox.minz) / 2.0f * finalzaxis;

		UnitBox ubox;

		for (auto ver : vertexset) {
			auto vh = mesh.vertex_handle(ver);
			auto pnt = mesh.point(vh);
			Eigen::Vector3d point{ pnt[0],pnt[1],pnt[2] };
			float xvalue = abs((point - newcentroid).dot(finalxaxis));
			float yvalue = abs((point - newcentroid).dot(finalyaxis));
			float zvalue = abs((point - newcentroid).dot(finalzaxis));
			ubox.refresh(xvalue, yvalue, zvalue);
		}
		float scalex = std::max(abs(ubox.minx), abs(ubox.maxx));
		float scaley = std::max(abs(ubox.miny), abs(ubox.maxy));
		float scalez = std::max(abs(ubox.minz), abs(ubox.maxz));
		float elemx = 1.0 / scalex;
		float elemy = 1.0 / scaley;
		float elemz = 1.0 / scalez;

		Eigen::Matrix3d rotationMatrix;
		rotationMatrix << finalxaxis, finalyaxis, finalzaxis;
		Eigen::Vector3d eulerXYZ = rotationMatrix.eulerAngles(0, 1, 2);

		SimpleCheck(eulerXYZ, finalxaxis, finalyaxis, finalzaxis);

		PackGEO& tempgeo = packexter[eid];
		tempgeo.needtransform = true;
		tempgeo.euler = eulerXYZ;
		tempgeo.translation = newcentroid;
		

		if (scalex < 1e-10) {
			scalex = 0.0f;
			tempgeo.isXPlatForm = true;
		}
		tempgeo.scaleelem = Eigen::Vector3d{ scalex,scaley,scalez };



		if (firstbox) {
			firstbox = false;
			MeshScaleBox = UnitBox{ scalex,scaley,scalez ,scalex,scaley,scalez };
			MeshTransBox = UnitBox{ float(newcentroid.x()),float(newcentroid.y()),float(newcentroid.z()),float(newcentroid.x()),float(newcentroid.y()),float(newcentroid.z()) };
		}
		else {
			if (MeshScaleBox.minx == 0.0f)
				MeshScaleBox.minx = MeshScaleBox.maxx;

			if (scalex != 0.0f)
				MeshScaleBox.refresh(scalex, scaley, scalez);
			else
				MeshScaleBox.refresh(MeshScaleBox.minx, scaley, scalez);

			MeshTransBox.refresh(float(newcentroid.x()), float(newcentroid.y()), float(newcentroid.z()));
		}

		double radx = eulerXYZ.x(), rady = eulerXYZ.y(), radz = eulerXYZ.z();
		tempgeo.rotatx = radx / PI * 255;
		tempgeo.rotaty = rady / PI * 255;
		tempgeo.rotatz = radz / PI * 255;

		for (auto& ver : vertexset) {
			auto vh = mesh.vertex_handle(ver);
			auto pnt = mesh.point(vh);
			Eigen::Vector3d point{ pnt[0],pnt[1],pnt[2] };
			//translation
			point = point - newcentroid;
			//rotation
			float valuex = point.dot(finalxaxis);
			float valuey = point.dot(finalyaxis);
			float valuez = point.dot(finalzaxis);
			//scale
			valuex *= elemx;
			valuey *= elemy;
			valuez *= elemz;
			assert(valuex <= 1.0);
			assert(valuey <= 1.0);
			assert(valuez <= 1.0);
			tempgeo.dataunit.emplace_back(valuex, valuey, valuez);
		}

		NextSimpleCheck(tempgeo, vertexset, mesh);

	}

	for (uint iid = 0; iid < targets.size(); ++iid) {
		auto vertexset = targets[iid].InternalLW.vertex;
		PackGEO& GEO = packinter[iid];
		if (vertexset.size() < 3) {
			GEO.needtransform = false;
			if(vertexset[0]!=EMPTYWIRE)
				for (auto& ver : vertexset)
				{
					auto vh = mesh.vertex_handle(ver);
					auto pnt = mesh.point(vh);
					GEO.nopackgeo.emplace_back(pnt[0], pnt[1], pnt[2]);
				}
			continue;
		}

		Eigen::MatrixXd WireVertices(3, vertexset.size());
		for (uint i = 0; i < vertexset.size(); ++i) {
			auto ver = vertexset[i];
			auto vh = mesh.vertex_handle(ver);
			auto pnt = mesh.point(vh);
			WireVertices.col(i) << pnt[0], pnt[1], pnt[2];
		}

		Eigen::Vector3d centroid = WireVertices.rowwise().mean().transpose();
		Eigen::MatrixXd centeredPoints = WireVertices.colwise() - centroid;
		// 执行SVD
		Eigen::JacobiSVD<Eigen::MatrixXd> svd(centeredPoints, Eigen::ComputeThinU | Eigen::ComputeThinV);
		// 获取平面的法向量（对应于最小奇异值的右奇异向量）
		Eigen::Vector3d normal = svd.matrixU().col(2).normalized();

		Eigen::Vector3d starty, startz;
		getStartAxis(starty, startz, normal);

		Box2d globox{ -1000000.0,-1000000.0,1000000.0,1000000.0 };
		uint glocnt = 0;
		float anglediff = 15.0f;

		for (uint anglecnt = 0; anglecnt < uint(90.0f / anglediff); anglecnt++) {
			float angle = anglecnt * anglediff;
			float radian = angle / 180.0f * PI;
			Eigen::Vector3d yaxis = cos(radian) * starty + sin(radian) * startz;
			Eigen::Vector3d zaxis = cos(radian) * startz - sin(radian) * starty;
			assert(yaxis.dot(zaxis) < 0.0001f);

			//计算2维包围盒
			Box2d box;

			for (auto ver : vertexset) {
				auto vh = mesh.vertex_handle(ver);
				auto pnt = mesh.point(vh);
				Eigen::Vector3d temppnt{ pnt[0],pnt[1],pnt[2] };
				Eigen::Vector3d p2p = temppnt - centroid;
				float yvalue = p2p.dot(yaxis);
				float zvalue = p2p.dot(zaxis);
				box.refresh(yvalue, zvalue);
			}
			if (box.getSumLength() < globox.getSumLength()) {
				globox = box;
				glocnt = anglecnt;
			}
		}

		Eigen::Vector3d finalxaxis = normal.normalized();
		Eigen::Vector3d finalyaxis = cos(glocnt * anglediff / 180.0 * PI) * starty + sin(glocnt * anglediff / 180.0 * PI) * startz;
		Eigen::Vector3d finalzaxis = cos(glocnt * anglediff / 180.0 * PI) * startz - sin(glocnt * anglediff / 180.0 * PI) * starty;
		Eigen::Vector3d newcentroid = centroid + (globox.maxy + globox.miny) / 2.0f * finalyaxis + (globox.maxz + globox.minz) / 2.0f * finalzaxis;

		UnitBox ubox;

		for (auto ver : vertexset) {
			auto vh = mesh.vertex_handle(ver);
			auto pnt = mesh.point(vh);
			Eigen::Vector3d point{ pnt[0],pnt[1],pnt[2] };
			float xvalue = abs((point - newcentroid).dot(finalxaxis));
			float yvalue = abs((point - newcentroid).dot(finalyaxis));
			float zvalue = abs((point - newcentroid).dot(finalzaxis));
			ubox.refresh(xvalue, yvalue, zvalue);
		}
		float scalex = std::max(abs(ubox.minx), abs(ubox.maxx));
		float scaley = std::max(abs(ubox.miny), abs(ubox.maxy));
		float scalez = std::max(abs(ubox.minz), abs(ubox.maxz));
		float elemx = 1.0 / scalex;
		float elemy = 1.0 / scaley;
		float elemz = 1.0 / scalez;

		Eigen::Matrix3d rotationMatrix;
		rotationMatrix << finalxaxis, finalyaxis, finalzaxis;
		Eigen::Vector3d eulerXYZ = rotationMatrix.eulerAngles(0, 1, 2);

		SimpleCheck(eulerXYZ, finalxaxis, finalyaxis, finalzaxis);

		PackGEO& tempgeo = packinter[iid];
		tempgeo.needtransform = true;
		tempgeo.euler = eulerXYZ;
		tempgeo.translation = newcentroid;

		if (scalex < 1e-10) {
			scalex = 0.0f;
			tempgeo.isXPlatForm = true;
		}
		tempgeo.scaleelem = Eigen::Vector3d{ scalex,scaley,scalez };
		if(scalex!=0.0f)
			MeshScaleBox.refresh(scalex, scaley, scalez);
		else
			MeshScaleBox.refresh(MeshScaleBox.minx, scaley, scalez);

		MeshTransBox.refresh(float(newcentroid.x()), float(newcentroid.y()), float(newcentroid.z()));

		double radx = eulerXYZ.x(), rady = eulerXYZ.y(), radz = eulerXYZ.z();
		tempgeo.rotatx = radx / PI * 255;
		tempgeo.rotaty = rady / PI * 255;
		tempgeo.rotatz = radz / PI * 255;


		for (auto& ver : vertexset) {
			auto vh = mesh.vertex_handle(ver);
			auto pnt = mesh.point(vh);
			Eigen::Vector3d point{ pnt[0],pnt[1],pnt[2] };
			//translation
			point = point - newcentroid;
			//rotation
			float valuex = point.dot(finalxaxis);
			float valuey = point.dot(finalyaxis);
			float valuez = point.dot(finalzaxis);
			//scale
			valuex *= elemx;
			valuey *= elemy;
			valuez *= elemz;
			assert(valuex <= 1.0);
			assert(valuey <= 1.0);
			assert(valuez <= 1.0);
			tempgeo.dataunit.emplace_back(valuex, valuey, valuez);
		}

		NextSimpleCheck(tempgeo, vertexset, mesh);
	}

	float xlengthtrans = MeshTransBox.maxx - MeshTransBox.minx;
	float ylengthtrans = MeshTransBox.maxy - MeshTransBox.miny;
	float zlengthtrans = MeshTransBox.maxz - MeshTransBox.minz;

	float undernumx = std::log(MeshScaleBox.maxx / MeshScaleBox.minx);
	float undernumy = std::log(MeshScaleBox.maxy / MeshScaleBox.miny);
	float undernumz = std::log(MeshScaleBox.maxz / MeshScaleBox.minz);


	for (uint eid = 0; eid < packexter.size(); eid++) {
		if (packexter[eid].needtransform == false)continue;
		auto& tempgeo = packexter[eid];
		tempgeo.translatex = static_cast<uchar>(255 * ((tempgeo.translation.x() - MeshTransBox.minx) / xlengthtrans));
		tempgeo.translatey = static_cast<uchar>(255 * ((tempgeo.translation.y() - MeshTransBox.miny) / ylengthtrans));
		tempgeo.translatez = static_cast<uchar>(255 * ((tempgeo.translation.z() - MeshTransBox.minz) / zlengthtrans));

		//tempgeo.scalex = ;
		if (tempgeo.isXPlatForm == false)
			tempgeo.scalex = static_cast<uchar>(255 * (std::log(tempgeo.scaleelem.x() / MeshScaleBox.minx) / undernumx));
		tempgeo.scaley = static_cast<uchar>(255 * (std::log(tempgeo.scaleelem.y() / MeshScaleBox.miny) / undernumy));
		tempgeo.scalez = static_cast<uchar>(255 * (std::log(tempgeo.scaleelem.z() / MeshScaleBox.minz) / undernumz));
	}

	for (uint iid = 0; iid < packinter.size(); iid++) {
		if (packinter[iid].needtransform == false)continue;
		auto& tempgeo = packinter[iid];
		tempgeo.translatex = static_cast<uchar>(255 * ((tempgeo.translation.x() - MeshTransBox.minx) / xlengthtrans));
		tempgeo.translatey = static_cast<uchar>(255 * ((tempgeo.translation.y() - MeshTransBox.miny) / ylengthtrans));
		tempgeo.translatez = static_cast<uchar>(255 * ((tempgeo.translation.z() - MeshTransBox.minz) / zlengthtrans));

		//tempgeo.scalex = ;
		if (tempgeo.isXPlatForm == false)
			tempgeo.scalex = static_cast<uchar>(255 * (std::log(tempgeo.scaleelem.x() / MeshScaleBox.minx) / undernumx));
		tempgeo.scaley = static_cast<uchar>(255 * (std::log(tempgeo.scaleelem.y() / MeshScaleBox.miny) / undernumy));
		tempgeo.scalez = static_cast<uchar>(255 * (std::log(tempgeo.scaleelem.z() / MeshScaleBox.minz) / undernumz));
	}
	
	//check by dequantize
	{
		for (uint eid = 0; eid < packexter.size(); eid++) {


			auto vertexset = gloEwires[eid].vertex;
			auto& tempgeo = packexter[eid];

			if (tempgeo.needtransform == false)
				continue;

			Eigen::Vector3d standardXAxis(1.0, 0.0, 0.0);
			Eigen::Vector3d standardYAxis(0.0, 1.0, 0.0);
			Eigen::Vector3d standardZAxis(0.0, 0.0, 1.0);

			Eigen::Matrix3d rotationMatrix;

			Eigen::Vector3d eulerxyz = Eigen::Vector3d{
				tempgeo.rotatx/255.0f*PI,tempgeo.rotaty/255.0f*PI,tempgeo.rotatz/255.0f*PI
			};
			rotationMatrix = Eigen::AngleAxisd(eulerxyz.x(), Eigen::Vector3d::UnitX())
				* Eigen::AngleAxisd(eulerxyz.y(), Eigen::Vector3d::UnitY())
				* Eigen::AngleAxisd(eulerxyz.z(), Eigen::Vector3d::UnitZ());

			Eigen::Vector3d newXAxis = rotationMatrix * standardXAxis;
			Eigen::Vector3d newYAxis = rotationMatrix * standardYAxis;
			Eigen::Vector3d newZAxis = rotationMatrix * standardZAxis;

			float descalex = MeshScaleBox.minx * std::pow(MeshScaleBox.maxx / MeshScaleBox.minx, tempgeo.scalex / 255.0f);
			float descaley = MeshScaleBox.miny * std::pow(MeshScaleBox.maxy / MeshScaleBox.miny, tempgeo.scaley / 255.0f);
			float descalez = MeshScaleBox.minz * std::pow(MeshScaleBox.maxz / MeshScaleBox.minz, tempgeo.scalez / 255.0f);

			float detransx = (tempgeo.translatex / 255.0f) * (MeshTransBox.maxx - MeshTransBox.minx) + MeshTransBox.minx;
			float detransy = (tempgeo.translatey / 255.0f) * (MeshTransBox.maxy - MeshTransBox.miny) + MeshTransBox.miny;
			float detransz = (tempgeo.translatez / 255.0f) * (MeshTransBox.maxz - MeshTransBox.minz) + MeshTransBox.minz;
			Eigen::Vector3d anotherdis{ detransx,detransy,detransz };

			for (uint i = 0; i < vertexset.size(); ++i) {
				auto pnt = mesh.point(mesh.vertex_handle(vertexset[i]));
				auto rightpnt = Eigen::Vector3d{ pnt[0],pnt[1],pnt[2] };

				Eigen::Vector3d rawdata = tempgeo.dataunit[i];
				Eigen::Vector3d temppnt = newXAxis * rawdata.x() * descalex + newYAxis * rawdata.y() * descaley + newZAxis * rawdata.z() * descalez;
				temppnt += anotherdis;
				//if((temppnt - rightpnt).norm()>0.1)
					std::cout<<"解压缩与实际值比较" << (temppnt - rightpnt).norm() << std::endl;
				
			}

			for (uint i = 0; i < vertexset.size(); ++i) {
				auto pnt = mesh.point(mesh.vertex_handle(vertexset[i]));
				auto rightpnt = Eigen::Vector3d{ pnt[0],pnt[1],pnt[2] };

				rightpnt -= anotherdis;
				float compressxvalue = rightpnt.dot(newXAxis);
				float compressyvalue = rightpnt.dot(newYAxis);
				float compresszvalue = rightpnt.dot(newZAxis);

				compressxvalue /= descalex;
				compressyvalue /= descaley;
				compresszvalue /= descalez;

				tempgeo.datatocompress.emplace_back(compressxvalue, compressyvalue, compresszvalue);
				//std::cout << (tempgeo.datatocompress[i] - tempgeo.dataunit[i]).length();

			}

		
		}

	}

	//check by dequantize
	//
	{
		for (uint iid = 0; iid < packinter.size(); iid++) {


			auto vertexset = targets[iid].InternalLW.vertex;
			auto& tempgeo = packinter[iid];

			if (tempgeo.needtransform == false)
				continue;

			Eigen::Vector3d standardXAxis(1.0, 0.0, 0.0);
			Eigen::Vector3d standardYAxis(0.0, 1.0, 0.0);
			Eigen::Vector3d standardZAxis(0.0, 0.0, 1.0);

			Eigen::Matrix3d rotationMatrix;

			Eigen::Vector3d eulerxyz = Eigen::Vector3d{
				tempgeo.rotatx / 255.0f * PI,tempgeo.rotaty / 255.0f * PI,tempgeo.rotatz / 255.0f * PI
			};
			rotationMatrix = Eigen::AngleAxisd(eulerxyz.x(), Eigen::Vector3d::UnitX())
				* Eigen::AngleAxisd(eulerxyz.y(), Eigen::Vector3d::UnitY())
				* Eigen::AngleAxisd(eulerxyz.z(), Eigen::Vector3d::UnitZ());

			Eigen::Vector3d newXAxis = rotationMatrix * standardXAxis;
			Eigen::Vector3d newYAxis = rotationMatrix * standardYAxis;
			Eigen::Vector3d newZAxis = rotationMatrix * standardZAxis;

			float descalex = MeshScaleBox.minx * std::pow(MeshScaleBox.maxx / MeshScaleBox.minx, tempgeo.scalex / 255.0f);
			float descaley = MeshScaleBox.miny * std::pow(MeshScaleBox.maxy / MeshScaleBox.miny, tempgeo.scaley / 255.0f);
			float descalez = MeshScaleBox.minz * std::pow(MeshScaleBox.maxz / MeshScaleBox.minz, tempgeo.scalez / 255.0f);

			float detransx = (tempgeo.translatex / 255.0f) * (MeshTransBox.maxx - MeshTransBox.minx) + MeshTransBox.minx;
			float detransy = (tempgeo.translatey / 255.0f) * (MeshTransBox.maxy - MeshTransBox.miny) + MeshTransBox.miny;
			float detransz = (tempgeo.translatez / 255.0f) * (MeshTransBox.maxz - MeshTransBox.minz) + MeshTransBox.minz;
			Eigen::Vector3d anotherdis{ detransx,detransy,detransz };

			for (uint i = 0; i < vertexset.size(); ++i) {
				auto pnt = mesh.point(mesh.vertex_handle(vertexset[i]));
				auto rightpnt = Eigen::Vector3d{ pnt[0],pnt[1],pnt[2] };

				Eigen::Vector3d rawdata = tempgeo.dataunit[i];
				Eigen::Vector3d temppnt = newXAxis * rawdata.x() * descalex + newYAxis * rawdata.y() * descaley + newZAxis * rawdata.z() * descalez;
				temppnt += anotherdis;
				//if((temppnt - rightpnt).norm()>0.1)
				std::cout << "解压缩与实际值比较" << (temppnt - rightpnt).norm() << std::endl;

			}


		}

	}

	//generate the presion value
	float error_threshold = 1e-6;
	for (uint eid = 0; eid < packexter.size(); eid++) {
		auto& tempgeo = packexter[eid];
		if (tempgeo.isXPlatForm == false)
			continue;


	}

}

void NewLWGenerator::fitPlane(const MyMesh& mesh, const std::vector<uint>& verset, Eigen::Vector3d& centroid, Eigen::Vector3d& normal)
{
	Eigen::MatrixXd WireVertices(3, verset.size());
	for (uint i = 0; i < verset.size(); ++i) {
		auto ver = verset[i];
		auto vh = mesh.vertex_handle(ver);
		auto pnt = mesh.point(vh);
		WireVertices.col(i) << pnt[0], pnt[1], pnt[2];
	}

	centroid = WireVertices.rowwise().mean().transpose();
	Eigen::MatrixXd centeredPoints = WireVertices.colwise() - centroid;
	Eigen::JacobiSVD<Eigen::MatrixXd> svd(centeredPoints, Eigen::ComputeThinU | Eigen::ComputeThinV);
	normal = svd.matrixU().col(2).normalized();
}

Box2d NewLWGenerator::getYZAxis(const MyMesh& mesh, const std::vector<uint>& verset, Eigen::Vector3d xaxis, Eigen::Vector3d& yaxis, Eigen::Vector3d& zaxis, float anglediff,Eigen::Vector3d centroid)
{
	Eigen::Vector3d starty, startz;
	getStartAxis(starty, startz, xaxis);

	std::vector<Eigen::Vector3d> pointset;
	for (const auto& ver : verset) {
		auto vh = mesh.vertex_handle(ver);
		auto pnt = mesh.point(vh);
		pointset.emplace_back(pnt[0] - centroid[0], pnt[1] - centroid[1], pnt[2] - centroid[2]);
	}


	Box2d globox{ -1000000.0,-1000000.0,1000000.0,1000000.0 };
	uint glocnt = 0;

	for (uint anglecnt = 0; anglecnt < uint(PI/2.0f / anglediff); anglecnt++) {
		float radian = anglecnt * anglediff;
		Eigen::Vector3d tempyaxis = cos(radian) * starty + sin(radian) * startz;
		Eigen::Vector3d tempzaxis = cos(radian) * startz - sin(radian) * starty;
		assert(tempyaxis.dot(tempzaxis) < 0.0001f);
		//计算2维包围盒
		Box2d box;

		for (auto pnt : pointset) {
			float yvalue = pnt.dot(tempyaxis);
			float zvalue = pnt.dot(tempzaxis);
			box.refresh(yvalue, zvalue);
		}
		if (box.getSumLength() < globox.getSumLength()) {
			globox = box;
			glocnt = anglecnt;
		}
	}

	float finalradian = glocnt * anglediff;
	yaxis = cos(finalradian) * starty + sin(finalradian) * startz;
	zaxis = cos(finalradian) * startz - sin(finalradian) * starty;

	return globox;
}

void NewLWGenerator::getScaleBox(const MyMesh& mesh, const std::vector<uint>& verset, Eigen::Vector3d centroid, UnitBox& scalebox, Eigen::Vector3d xaxis, Eigen::Vector3d& yaxis, Eigen::Vector3d& zaxis)
{
	for (const auto& ver : verset) {
		auto vh = mesh.vertex_handle(ver);
		auto pnt = mesh.point(vh);
		Eigen::Vector3d temppoint{ pnt[0] - centroid[0],pnt[1] - centroid[1],pnt[2] - centroid[2] };
		float xvalue = temppoint.dot(xaxis);
		float yvalue = temppoint.dot(yaxis);
		float zvalue = temppoint.dot(zaxis);

		scalebox.refresh(std::abs(xvalue), std::abs(yvalue), std::abs(zvalue));
	}
}

void NewLWGenerator::NewVertexQuantization(const MyMesh& mesh)
{
	packexter.resize(gloEwires.size());
	packinter.resize(targets.size());


	bool firstadd = true;

	for (uint eid = 0; eid < packexter.size(); ++eid) {
		//个数小于3的就不处理
		if (gloEwires[eid].vertex.size() < 3) {
			PackNoPlane(packexter[eid], mesh, gloEwires[eid]);
			continue;
		}

		auto& tempgeo = packexter[eid];
		auto& ewire = gloEwires[eid];

		Eigen::Vector3d centroid, normal,newcentroid;
		Eigen::Vector3d finalyaxis, finalzaxis;
		Box2d globox;
		UnitBox scalebox;

		fitPlane(mesh, ewire.vertex, centroid, normal);
		globox = getYZAxis(mesh, ewire.vertex, normal,finalyaxis,finalzaxis,15.0/180.0*PI,centroid);
		newcentroid = centroid + (globox.maxy + globox.miny) / 2.0f * finalyaxis + (globox.maxz + globox.minz) / 2.0f * finalzaxis;
		getScaleBox(mesh, ewire.vertex, newcentroid, scalebox,normal,finalyaxis,finalzaxis);
		
		//pack tempgeo
		tempgeo.needtransform = true;
		tempgeo.translation = newcentroid;
		tempgeo.scaleelem = Eigen::Vector3d{ scalebox.maxx,scalebox.maxy,scalebox.maxz };
		Eigen::Matrix3d rotationMatrix;
		rotationMatrix << normal, finalyaxis, finalzaxis;
		Eigen::Vector3d eulerXYZ = rotationMatrix.eulerAngles(0, 1, 2);
		tempgeo.euler = eulerXYZ;

		//Eigen::AngleAxisd rotationX(eulerXYZ[0], Eigen::Vector3d::UnitX());
		//Eigen::AngleAxisd rotationY(eulerXYZ[1], Eigen::Vector3d::UnitY());
		//Eigen::AngleAxisd rotationZ(eulerXYZ[2], Eigen::Vector3d::UnitZ());
		//Eigen::Matrix3d newrotationMatrix = rotationX.matrix() * rotationY.matrix() * rotationZ.matrix();
		//相当于一个reverse rotation matrix

		//	x - axis
		//  y - axis  --> transpose --> newratationmatrix (x*y*z)
		//  z - axis

		//refresh meshtranbox ,meshscalebox
		NewSimpleCheck(tempgeo, mesh, ewire.vertex);

		if (firstadd) {
			MeshTransBox = UnitBox{ float(newcentroid[0]),float(newcentroid[1]),float(newcentroid[2])
				,float(newcentroid[0]),float(newcentroid[1]),float(newcentroid[2]) };
			MeshScaleBox = UnitBox{
				scalebox.maxx,scalebox.maxy,scalebox.maxz,
				scalebox.maxx,scalebox.maxy,scalebox.maxz };
			firstadd = false;
		}
		else {
			MeshTransBox.refresh(newcentroid[0], newcentroid[1], newcentroid[2]);
			if (scalebox.maxx < 1e-10)
				scalebox.maxx = MeshScaleBox.minx;
			MeshScaleBox.refresh(scalebox.maxx, scalebox.maxy, scalebox.maxz);
		}
	}

	for (uint iid = 0; iid < packinter.size(); ++iid) {
		//个数小于3的就不处理
		if (targets[iid].InternalLW.vertex.size() < 3) {
			PackNoPlane(packinter[iid], mesh, targets[iid].InternalLW);
			continue;
		}

		auto& tempgeo = packinter[iid];
		auto& ewire = targets[iid].InternalLW;

		Eigen::Vector3d centroid, normal, newcentroid;
		Eigen::Vector3d finalyaxis, finalzaxis;
		Box2d globox;
		UnitBox scalebox;

		fitPlane(mesh, ewire.vertex, centroid, normal);
		globox = getYZAxis(mesh, ewire.vertex, normal, finalyaxis, finalzaxis, 15.0 / 180.0 * PI, centroid);
		newcentroid = centroid + (globox.maxy + globox.miny) / 2.0f * finalyaxis + (globox.maxz + globox.minz) / 2.0f * finalzaxis;
		getScaleBox(mesh, ewire.vertex, newcentroid, scalebox, normal, finalyaxis, finalzaxis);

		//pack tempgeo
		tempgeo.needtransform = true;
		tempgeo.translation = newcentroid;
		tempgeo.scaleelem = Eigen::Vector3d{ scalebox.maxx,scalebox.maxy,scalebox.maxz };
		Eigen::Matrix3d rotationMatrix;
		rotationMatrix << normal, finalyaxis, finalzaxis;
		Eigen::Vector3d eulerXYZ = rotationMatrix.eulerAngles(0, 1, 2);
		tempgeo.euler = eulerXYZ;

		//refresh meshtranbox ,meshscalebox
		NewSimpleCheck(tempgeo, mesh, ewire.vertex);

		MeshTransBox.refresh(newcentroid[0], newcentroid[1], newcentroid[2]);
		if (scalebox.maxx < 1e-10)
			scalebox.maxx = MeshScaleBox.minx;
		MeshScaleBox.refresh(scalebox.maxx, scalebox.maxy, scalebox.maxz);
	}
	UcharRadGen();

	//for (uint eid = 0; eid < gloEwires.size(); ++eid)
	//	if(packexter[eid].needtransform == true)
	//		NewSimpleCheck(packexter[eid], mesh, gloEwires[eid].vertex);
	TranslatePack(mesh, MeshTransBox);

	FinalScaleGen(mesh);
	//pack data to compress;
	CheckByDequantize(mesh);
	//vertexbitGen();

	return;
}

void NewLWGenerator::PackNoPlane(PackGEO& tempgeo, const MyMesh& mesh, const MeshletBody& ewire)
{
	tempgeo.needtransform = false;
	for (const auto& ver : ewire.vertex) {
		if (ver == EMPTYWIRE)
			return;
		auto pnt = mesh.point(mesh.vertex_handle(ver));
		tempgeo.nopackgeo.emplace_back(pnt[0], pnt[1], pnt[2]);
	}
}

void NewLWGenerator::PackGPULW(const MyMesh& mesh)
{
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<float> dis(0.0, 1.0);

	std::vector<MeshletDes> records;
	records.resize(targets.size());

	//pack inter geo
	for (uint mid = 0; mid < targets.size(); ++mid) {
		auto& intermesh = targets[mid].InternalLW;
		uint vertexnum;
		if (intermesh.vertex[0] == EMPTYWIRE)
			vertexnum = 0;
		else
			vertexnum = intermesh.vertex.size();

		records[mid].ewirenum = static_cast<uchar>(targets[mid].ewireidx.size());
		records[mid].irrnum = static_cast<uchar>(intermesh.irregular.size()/3);
		records[mid].numvertex = static_cast<uchar>(internalbuilders[mid].vertices.size());

		records[mid].color[0] = static_cast<uchar>(dis(gen) * 255);
		records[mid].color[1] = static_cast<uchar>(dis(gen) * 255);
		records[mid].color[2] = static_cast<uchar>(dis(gen) * 255);

		records[mid].numinver = vertexnum;
		//records[mid].useless2 = 0;


		uint location = intergeo.size();
		records[mid].ingeolocation = location;
		records[mid].inconlocation = intercon.size()/4;
		
		if (vertexnum != 0) {
			for (auto& ver : intermesh.vertex) {
				auto vh = mesh.vertex_handle(ver);
				auto pnt = mesh.point(vh);
				intergeo.push_back(pnt[0]);
				intergeo.push_back(pnt[1]);
				intergeo.push_back(pnt[2]);
			}

			assert(intermesh.left.size() == intermesh.vertex.size());
			assert(intermesh.left.size() == intermesh.right.size());

			for (auto& elem : intermesh.left)
				intercon.push_back(elem);
			for (auto& elem : intermesh.right)
				intercon.push_back(elem);
		}
		
		for (auto& elem : intermesh.irregular)
			intercon.push_back(elem);

		while (intercon.size()%4!=0)
		{
			intercon.push_back(0);
		}
	}

	assert(intergeo.size() < (0xFFFFFFFF));
	assert(intercon.size() < 0xFFFFFFFF);

	//pack globol ewire
	std::vector<uint> Econloc;
	std::vector<uint> Egeoloc;
	for (auto& wire : gloEwires) {
		Egeoloc.push_back(extergeo.size());
		Econloc.push_back(extercon.size()/4);

		for (auto& ver : wire.vertex) {
			auto vh = mesh.vertex_handle(ver);
			auto pnt = mesh.point(vh);
			extergeo.push_back(pnt[0]);
			extergeo.push_back(pnt[1]);
			extergeo.push_back(pnt[2]);
		}

		for (auto& elem : wire.left)
			extercon.push_back(elem);
		for (uint i = 0; i < wire.right.size(); ++i)
			extercon.push_back(wire.right[wire.right.size() - i - 1]);

		while (extercon.size() % 4 != 0)
		{
			extercon.push_back(0);
		}

		assert(wire.right.size() == wire.left.size());
		assert(wire.right.size() == wire.vertex.size());
	}

	assert(extergeo.size() < (0xFFFFFFFF));
	assert(extercon.size() < 0xFFFFFFFF);

	for (uint mid = 0; mid < targets.size(); ++mid) {
		int cnt = 0;
		for (auto& id : targets[mid].ewireidx) {

			uint location = Egeoloc[id];
			uchar num = static_cast<uchar>(gloEwires[id].vertex.size());

			if (targets[mid].reverse[cnt] == true)
				num = num | 0x80;
			records[mid].numexver.push_back(num);
			records[mid].exgeolocation.push_back(location);
			records[mid].exconlocation.push_back(Econloc[id]);
			cnt++;
		}
	}

	for (uint mid = 0; mid < targets.size(); ++mid) {
		DesLoc.push_back(Desinfo.size());

		while ((records[mid].numexver.size() % 4) != 0)
			records[mid].numexver.push_back(0);
		records[mid].ingeostart = static_cast<uchar>(2 + records[mid].numexver.size() / 4);

		Desinfo.push_back(PackChar4Uint(records[mid].ewirenum, records[mid].color[0],
			records[mid].color[1], records[mid].color[2]));

		Desinfo.push_back(PackChar4Uint(records[mid].irrnum, records[mid].numvertex,
			records[mid].numinver, records[mid].ingeostart));

		for (int i = 0; i < records[mid].numexver.size(); i += 4) {
			Desinfo.push_back(PackChar4Uint(records[mid].numexver[i], records[mid].numexver[i+1],
				records[mid].numexver[i+2], records[mid].numexver[i+3]));
		}

		Desinfo.push_back(records[mid].ingeolocation);
		Desinfo.push_back(records[mid].inconlocation);
		for (auto& elem : records[mid].exgeolocation)
			Desinfo.push_back(elem);
		for (auto& elem : records[mid].exconlocation)
			Desinfo.push_back(elem);
	}

	assert(Desinfo.size() < 0xFFFFFFFF);

	assert(extercon.size() % 4 == 0);
	assert(intercon.size() % 4 == 0);
	for (uint i = 0; i < extercon.size(); i += 4)
		newextercon.push_back(PackChar4Uint(extercon[i], extercon[i + 1], extercon[i + 2], extercon[i + 3]));

	for (uint i = 0; i < intercon.size(); i += 4)
		newintercon.push_back(PackChar4Uint(intercon[i], intercon[i + 1], intercon[i + 2], intercon[i + 3]));

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
void NewLWGenerator::getStartAxis(Eigen::Vector3d& yaxis, Eigen::Vector3d& zaxis, const Eigen::Vector3d& normal)
{
	Eigen::Vector3d temp{ 1.0f,0.0f,0.0f };
	if (temp.dot(normal) == 1.0f)
		temp = Eigen::Vector3d{ 0.0f,1.0f,0.0f };

	yaxis = temp.cross(normal).normalized();
	zaxis = normal.cross(yaxis).normalized();

}
bool NewLWGenerator::SimpleCheck(Eigen::Vector3d eulerxyz, Eigen::Vector3d newx, Eigen::Vector3d newy, Eigen::Vector3d newz)
{
	
	Eigen::Vector3d standardXAxis(1.0, 0.0, 0.0);
	Eigen::Vector3d standardYAxis(0.0, 1.0, 0.0);
	Eigen::Vector3d standardZAxis(0.0, 0.0, 1.0);

	Eigen::Matrix3d rotationMatrix;
	rotationMatrix = Eigen::AngleAxisd(eulerxyz.x(), Eigen::Vector3d::UnitX())
		* Eigen::AngleAxisd(eulerxyz.y(), Eigen::Vector3d::UnitY())
		* Eigen::AngleAxisd(eulerxyz.z(), Eigen::Vector3d::UnitZ());

	Eigen::Vector3d newXAxis = rotationMatrix * standardXAxis;
	Eigen::Vector3d newYAxis = rotationMatrix * standardYAxis;
	Eigen::Vector3d newZAxis = rotationMatrix * standardZAxis;

	std::cout<<eulerxyz.x()/PI*180<< ' ' << (newx - newXAxis).norm() << std::endl;
	std::cout << eulerxyz.y() / PI * 180 << ' ' << (newy - newYAxis).norm() << std::endl;
	std::cout << eulerxyz.z() / PI * 180 << ' ' << (newz - newZAxis).norm() << std::endl;

	return true;
}
void NewLWGenerator::NextSimpleCheck(const PackGEO& tempgeo,const std::vector<uint>& vertexvec, const MyMesh& mesh)
{
	Eigen::Vector3d standardXAxis(1.0, 0.0, 0.0);
	Eigen::Vector3d standardYAxis(0.0, 1.0, 0.0);
	Eigen::Vector3d standardZAxis(0.0, 0.0, 1.0);

	Eigen::Matrix3d rotationMatrix;

	Eigen::Vector3d eulerxyz = tempgeo.euler;
	rotationMatrix = Eigen::AngleAxisd(eulerxyz.x(), Eigen::Vector3d::UnitX())
		* Eigen::AngleAxisd(eulerxyz.y(), Eigen::Vector3d::UnitY())
		* Eigen::AngleAxisd(eulerxyz.z(), Eigen::Vector3d::UnitZ());

	Eigen::Vector3d newXAxis = rotationMatrix * standardXAxis;
	Eigen::Vector3d newYAxis = rotationMatrix * standardYAxis;
	Eigen::Vector3d newZAxis = rotationMatrix * standardZAxis;

	for (uint i = 0; i < vertexvec.size(); ++i) {
		auto vh = mesh.vertex_handle(vertexvec[i]);
		auto pnt = mesh.point(vh);
		Eigen::Vector3d rightpoint{ pnt[0],pnt[1],pnt[2] };

		auto uintdata = tempgeo.dataunit[i];
		Eigen::Vector3d cmppoint{
			uintdata[0] * tempgeo.scaleelem.x(),
			uintdata[1] * tempgeo.scaleelem.y(),
			uintdata[2] * tempgeo.scaleelem.z(),
		};
		cmppoint = cmppoint[0] * newXAxis + cmppoint[1] * newYAxis+ cmppoint[2]* newZAxis+tempgeo.translation;
		auto distance = (cmppoint - rightpoint).norm();
		//std::cout<<"check in 恢复data" << distance  << std::endl;
		assert(distance < 1e-8);
	}
}

void NewLWGenerator::NewSimpleCheck(PackGEO& tempgeo, const MyMesh& mesh, const std::vector<uint>& verset)
{
	Eigen::Vector3d eulerXYZ = tempgeo.euler;
	Eigen::AngleAxisd rotationX(eulerXYZ[0], Eigen::Vector3d::UnitX());
	Eigen::AngleAxisd rotationY(eulerXYZ[1], Eigen::Vector3d::UnitY());
	Eigen::AngleAxisd rotationZ(eulerXYZ[2], Eigen::Vector3d::UnitZ());
	//dequantize
	Eigen::Matrix3d rotationMatrix = rotationX.matrix() * rotationY.matrix() * rotationZ.matrix();
	Eigen::DiagonalMatrix<double, 3> scaleMatirx;
	scaleMatirx.diagonal() << tempgeo.scaleelem[0], tempgeo.scaleelem[1], tempgeo.scaleelem[2];

	//quanitize
	Eigen::Matrix3d revrotat = rotationMatrix.transpose();
	Eigen::DiagonalMatrix<double, 3> revscale;
	revscale.diagonal() << 1.0/tempgeo.scaleelem[0], 1.0/tempgeo.scaleelem[1], 1.0/tempgeo.scaleelem[2];


	Eigen::Matrix3d com = revscale * revrotat ;
	Eigen::Matrix3d decom =   rotationMatrix* scaleMatirx;
	//
	for (const auto& ver : verset) {
		auto vh = mesh.vertex_handle(ver);
		auto pnt = mesh.point(vh);
		auto rawdata = Eigen::Vector3d{ pnt[0],pnt[1],pnt[2] };

		Eigen::Vector3d unitdata = com * (rawdata - tempgeo.translation);
		tempgeo.dataunit.push_back(unitdata);
		Eigen::Vector3d cmpdata = decom * unitdata + tempgeo.translation;

		double diff = (cmpdata - rawdata).norm();
		assert(diff < 1e-12);
		std::cout << "distance cmp and raw data" << diff << std::endl;
	}

}

void NewLWGenerator::CheckByDequantize(const MyMesh& mesh)
{
	for (uint eid = 0; eid < packexter.size(); ++eid) {
		auto& tempgeo = packexter[eid];
		if (tempgeo.needtransform == false)continue;
		tempgeo.datatocompress.clear();

		Eigen::Vector3d eulerXYZ{
			float(tempgeo.rotatx) / 255.0 * TWOPI,
			float(tempgeo.rotaty) / 255.0 * TWOPI,
			float(tempgeo.rotatz) / 255.0 * TWOPI
		};
		Eigen::AngleAxisd rotationX(eulerXYZ[0], Eigen::Vector3d::UnitX());
		Eigen::AngleAxisd rotationY(eulerXYZ[1], Eigen::Vector3d::UnitY());
		Eigen::AngleAxisd rotationZ(eulerXYZ[2], Eigen::Vector3d::UnitZ());
		//quantize
		Eigen::Matrix3d rotationMatrix = (rotationX.matrix() * rotationY.matrix() * rotationZ.matrix()).transpose();
		Eigen::Vector3d tranvec = {
			MeshTransBox.minx + (tempgeo.translatex / 255.0) * (MeshTransBox.maxx - MeshTransBox.minx),
			MeshTransBox.miny + (tempgeo.translatey / 255.0) * (MeshTransBox.maxy - MeshTransBox.miny),
			MeshTransBox.minz + (tempgeo.translatez / 255.0) * (MeshTransBox.maxz - MeshTransBox.minz)
		};
		float descalex = MeshScaleBox.minx * std::pow(MeshScaleBox.maxx / MeshScaleBox.minx, tempgeo.scalex / 255.0f);
		float descaley = MeshScaleBox.miny * std::pow(MeshScaleBox.maxy / MeshScaleBox.miny, tempgeo.scaley / 255.0f);
		float descalez = MeshScaleBox.minz * std::pow(MeshScaleBox.maxz / MeshScaleBox.minz, tempgeo.scalez / 255.0f);

		Eigen::DiagonalMatrix<double, 3> revscale,scalema;
		revscale.diagonal() << 1.0 / descalex, 1.0 / descaley, 1.0 / descalez;
		scalema.diagonal() << descalex, descaley, descalez;
		Eigen::MatrixXd compressmatrix = revscale * rotationMatrix;
		Eigen::MatrixXd dequanmatrix = rotationMatrix.transpose() * scalema;

		auto vertexset = gloEwires[eid].vertex;
		for (auto ver : vertexset) {
			auto pnt = mesh.point(mesh.vertex_handle(ver));
			auto rawdata = Eigen::Vector3d{ pnt[0],pnt[1],pnt[2] };

			Eigen::Vector3d result = compressmatrix * (rawdata - tranvec);
			tempgeo.datatocompress.push_back(result);
			Eigen::Vector3d dequandata = dequanmatrix * result + tranvec;
			std::cout<<"差距是 " << (dequandata - rawdata).norm() << std::endl;
		}


	}

}

void NewLWGenerator::FinalScaleGen(const MyMesh& mesh)
{
	std::vector<float> xvalues;
	std::vector<float> yvalues;
	std::vector<float> zvalues;

	for (uint eid = 0; eid < packexter.size(); ++eid) {
		auto& tempgeo = packexter[eid];
		if (tempgeo.needtransform == false)continue;

		Eigen::Vector3d eulerXYZ{
			float(tempgeo.rotatx) / 255.0 * TWOPI,
			float(tempgeo.rotaty) / 255.0 * TWOPI,
			float(tempgeo.rotatz) / 255.0 * TWOPI
		};
		Eigen::AngleAxisd rotationX(eulerXYZ[0], Eigen::Vector3d::UnitX());
		Eigen::AngleAxisd rotationY(eulerXYZ[1], Eigen::Vector3d::UnitY());
		Eigen::AngleAxisd rotationZ(eulerXYZ[2], Eigen::Vector3d::UnitZ());
		//quantize
		Eigen::Matrix3d rotationMatrix = (rotationX.matrix() * rotationY.matrix() * rotationZ.matrix()).transpose();
		Eigen::Vector3d tranvec = {
			MeshTransBox.minx + (tempgeo.translatex / 255.0) * (MeshTransBox.maxx - MeshTransBox.minx),
			MeshTransBox.miny + (tempgeo.translatey / 255.0) * (MeshTransBox.maxy - MeshTransBox.miny),
			MeshTransBox.minz + (tempgeo.translatez / 255.0) * (MeshTransBox.maxz - MeshTransBox.minz)
		};

		auto vertexset = gloEwires[eid].vertex;
		float tempx = 0, tempy = 0,tempz = 0;
		for (auto ver : vertexset) {
			auto pnt = mesh.point(mesh.vertex_handle(ver));
			Eigen::Vector3d point{ pnt[0],pnt[1],pnt[2] };
			point -= tranvec;
			Eigen::Vector3d scalevec = rotationMatrix * point;
			if (scalevec.x() > tempx)tempx = scalevec.x();
			if (scalevec.y() > tempy)tempy = scalevec.y();
			if (scalevec.z() > tempz)tempz = scalevec.z();
		}
		xvalues.push_back(tempx); yvalues.push_back(tempy); zvalues.push_back(tempz);
	}

	for (uint iid = 0; iid < packinter.size(); ++iid) {
		auto& tempgeo = packinter[iid];
		if (tempgeo.needtransform == false)continue;

		Eigen::Vector3d eulerXYZ{
			float(tempgeo.rotatx) / 255.0 * TWOPI,
			float(tempgeo.rotaty) / 255.0 * TWOPI,
			float(tempgeo.rotatz) / 255.0 * TWOPI
		};
		Eigen::AngleAxisd rotationX(eulerXYZ[0], Eigen::Vector3d::UnitX());
		Eigen::AngleAxisd rotationY(eulerXYZ[1], Eigen::Vector3d::UnitY());
		Eigen::AngleAxisd rotationZ(eulerXYZ[2], Eigen::Vector3d::UnitZ());
		//quantize
		Eigen::Matrix3d rotationMatrix = (rotationX.matrix() * rotationY.matrix() * rotationZ.matrix()).transpose();
		Eigen::Vector3d tranvec = {
			MeshTransBox.minx + (tempgeo.translatex / 255.0) * (MeshTransBox.maxx - MeshTransBox.minx),
			MeshTransBox.miny + (tempgeo.translatey / 255.0) * (MeshTransBox.maxy - MeshTransBox.miny),
			MeshTransBox.minz + (tempgeo.translatez / 255.0) * (MeshTransBox.maxz - MeshTransBox.minz)
		};

		auto vertexset = targets[iid].InternalLW.vertex;
		float tempx = 0, tempy = 0, tempz = 0;
		for (auto ver : vertexset) {
			auto pnt = mesh.point(mesh.vertex_handle(ver));
			Eigen::Vector3d point{ pnt[0],pnt[1],pnt[2] };
			point -= tranvec;
			Eigen::Vector3d scalevec = rotationMatrix * point;
			if (scalevec.x() > tempx)tempx = scalevec.x();
			if (scalevec.y() > tempy)tempy = scalevec.y();
			if (scalevec.z() > tempz)tempz = scalevec.z();
		}
		xvalues.push_back(tempx); yvalues.push_back(tempy); zvalues.push_back(tempz);
	}

	MeshScaleBox = UnitBox{ xvalues.front(),yvalues.front(),zvalues.front(),xvalues.front(),yvalues.front(),zvalues.front() };
	//pack 部分
	for (auto elem : xvalues) {
		if ((elem > 1e-10) && elem < MeshScaleBox.minx)MeshScaleBox.minx = elem;
		if (elem > MeshScaleBox.maxx)MeshScaleBox.maxx = elem;
	}
	for (auto elem : yvalues) {
		if ((elem > 1e-10) && elem < MeshScaleBox.miny)MeshScaleBox.miny = elem;
		if (elem > MeshScaleBox.maxy)MeshScaleBox.maxy = elem;
	}
	for (auto elem : zvalues) {
		if ((elem > 1e-10) && elem < MeshScaleBox.minz)MeshScaleBox.minz = elem;
		if (elem > MeshScaleBox.maxz)MeshScaleBox.maxz = elem;
	}

	uint meshletid = 0;

	float undernumx = std::log(MeshScaleBox.maxx / MeshScaleBox.minx);
	float undernumy = std::log(MeshScaleBox.maxy / MeshScaleBox.miny);
	float undernumz = std::log(MeshScaleBox.maxz / MeshScaleBox.minz);

	for (uint eid = 0; eid < packexter.size(); ++eid) {
		auto& tempgeo = packexter[eid];
		if (tempgeo.needtransform == false)continue;

		if (xvalues[meshletid] < 1e-10)tempgeo.isXPlatForm = true;
		else tempgeo.scalex = static_cast<uchar>(255 * (std::log(xvalues[meshletid] / MeshScaleBox.minx) / undernumx));
		tempgeo.scaley = static_cast<uchar>(255 * (std::log(yvalues[meshletid] / MeshScaleBox.miny) / undernumy));
		tempgeo.scalez = static_cast<uchar>(255 * (std::log(zvalues[meshletid] / MeshScaleBox.minz) / undernumz));
		meshletid++;
	}
	for (uint iid = 0; iid < packinter.size(); ++iid) {
		auto& tempgeo = packinter[iid];
		if (tempgeo.needtransform == false)continue;
		if (xvalues[meshletid] < 1e-10)tempgeo.isXPlatForm = true;
		else tempgeo.scalex = static_cast<uchar>(255 * (std::log(xvalues[meshletid] / MeshScaleBox.minx) / undernumx));
		tempgeo.scaley = static_cast<uchar>(255 * (std::log(yvalues[meshletid] / MeshScaleBox.miny) / undernumy));
		tempgeo.scalez = static_cast<uchar>(255 * (std::log(zvalues[meshletid] / MeshScaleBox.minz) / undernumz));
		meshletid++;
	}

}

void NewLWGenerator::UcharRadGen()
{
	for (auto& tempgeo : packexter) {
		if (tempgeo.needtransform == false)
			continue;
		auto& euler = tempgeo.euler;
		if (euler.x() < 0)euler.x() += 2 * PI;
		if (euler.y() < 0)euler.y() += 2 * PI;
		if (euler.z() < 0)euler.z() += 2 * PI;

		assert(2 * PI >= euler.x());
		assert(2 * PI >= euler.y());
		assert(2 * PI >= euler.z());

		float value;
		value = 255.0 * euler.x() / (PI * 2);
		tempgeo.rotatx = myclamp(value, 0, 255);

		value = 255.0 * euler.y() / (PI * 2);
		tempgeo.rotaty = myclamp(value, 0, 255);

		value = 255.0 * euler.z() / (PI * 2);
		tempgeo.rotatz = myclamp(value, 0, 255);
	}

	for (auto& tempgeo : packinter) {
		if (tempgeo.needtransform == false)
			continue;
		auto& euler = tempgeo.euler;
		if (euler.x() < 0)euler.x() += 2 * PI;
		if (euler.y() < 0)euler.y() += 2 * PI;
		if (euler.z() < 0)euler.z() += 2 * PI;

		assert(2*PI >= euler.x());
		assert(2*PI >= euler.y());
		assert(2*PI >= euler.z());

		float value;
		value = 255.0 * euler.x() / (PI*2);
		tempgeo.rotatx = myclamp(value, 0, 255);

		value = 255.0 * euler.y() / (PI*2);
		tempgeo.rotaty = myclamp(value, 0, 255);

		value = 255.0 * euler.z() / (PI*2);
		tempgeo.rotatz = myclamp(value, 0, 255);

	}
}

void NewLWGenerator::TranslatePack(const MyMesh& mesh, UnitBox TranslateBox)
{

	//更新每一个meshlet的diff

	for (uint eid = 0; eid < packexter.size(); ++eid) {
		auto& tempgeo = packexter[eid];
		if (tempgeo.needtransform == false)continue;
		//压缩角度
		Eigen::Vector3d eulerXYZ{
			float(tempgeo.rotatx) / 255.0 * TWOPI,
			float(tempgeo.rotaty) / 255.0 * TWOPI,
			float(tempgeo.rotatz) / 255.0 * TWOPI
		};
		Eigen::AngleAxisd rotationX(eulerXYZ[0], Eigen::Vector3d::UnitX());
		Eigen::AngleAxisd rotationY(eulerXYZ[1], Eigen::Vector3d::UnitY());
		Eigen::AngleAxisd rotationZ(eulerXYZ[2], Eigen::Vector3d::UnitZ());
		//dequantize
		Eigen::Matrix3d rotationMatrix = rotationX.matrix() * rotationY.matrix() * rotationZ.matrix();
		Eigen::DiagonalMatrix<double, 3> scaleMatirx;
		scaleMatirx.diagonal() << tempgeo.scaleelem[0], tempgeo.scaleelem[1], tempgeo.scaleelem[2];
		Eigen::Matrix3d decom = rotationMatrix * scaleMatirx;

		//
		float totalabs = 0;
		Eigen::Vector3d diffmean{ 0,0,0 };

		for (uint i = 0; i < gloEwires[eid].vertex.size(); ++i) {
			auto pnt = mesh.point(mesh.vertex_handle(gloEwires[eid].vertex[i]));
			Eigen::Vector3d rawdata{ pnt[0],pnt[1],pnt[2] };
			Eigen::Vector3d cmpdata = decom * tempgeo.dataunit[i] + tempgeo.translation;
			//std::cout << "distance after compact rad " << (cmpdata - rawdata).norm() << std::endl;
			totalabs += std::abs((rawdata - cmpdata).norm());
			diffmean += rawdata - cmpdata;
		}
		std::cout << "更新前 diff " << totalabs << std::endl;
		diffmean /= gloEwires[eid].vertex.size();

		Eigen::Vector3d newtranslation = diffmean + tempgeo.translation;
		totalabs = 0;

		for (uint i = 0; i < gloEwires[eid].vertex.size(); ++i) {
			auto pnt = mesh.point(mesh.vertex_handle(gloEwires[eid].vertex[i]));
			Eigen::Vector3d rawdata{ pnt[0],pnt[1],pnt[2] };
			Eigen::Vector3d cmpdata = decom * tempgeo.dataunit[i] + newtranslation;
			//std::cout << "distance after compact rad " << (cmpdata - rawdata).norm() << std::endl;
			totalabs += std::abs((rawdata - cmpdata).norm());
		}
		std::cout << "更新后 diff " << totalabs << std::endl;
		//绝大多数的error变少了但是也有少部分error变多了
	}

	for (uint eid = 0; eid < packexter.size(); ++eid) {
		auto& tempgeo = packexter[eid];
		if (tempgeo.needtransform == false)continue;
		float tempx = (tempgeo.translation.x() - TranslateBox.minx) / (TranslateBox.maxx - TranslateBox.minx);
		float tempy = (tempgeo.translation.y() - TranslateBox.miny) / (TranslateBox.maxy - TranslateBox.miny);
		float tempz = (tempgeo.translation.z() - TranslateBox.minz) / (TranslateBox.maxz - TranslateBox.minz);
		tempgeo.translatex = myclamp(tempx * 255, 0, 255);
		tempgeo.translatey = myclamp(tempy * 255, 0, 255);
		tempgeo.translatez = myclamp(tempz * 255, 0, 255);
	}

	for (uint iid = 0; iid < packinter.size(); ++iid) {
		auto& tempgeo = packinter[iid];
		if (tempgeo.needtransform == false)continue;
		float tempx = (tempgeo.translation.x() - TranslateBox.minx) / (TranslateBox.maxx - TranslateBox.minx);
		float tempy = (tempgeo.translation.y() - TranslateBox.miny) / (TranslateBox.maxy - TranslateBox.miny);
		float tempz = (tempgeo.translation.z() - TranslateBox.minz) / (TranslateBox.maxz - TranslateBox.minz);
		tempgeo.translatex = myclamp(tempx * 255, 0, 255);
		tempgeo.translatey = myclamp(tempy * 255, 0, 255);
		tempgeo.translatez = myclamp(tempz * 255, 0, 255);
	}

}

uchar NewLWGenerator::myclamp(float value, uchar lowbound, uchar highbound)
{
	if (value > highbound)return highbound;
	else if (value < lowbound)return lowbound;
	return uchar(uint(value));
}

void NewLWGenerator::ExportFile(const std::string& filename)
{
	size_t lastDotpos = filename.find_last_of('.');
	std::string partfilename = filename.substr(0, lastDotpos + 1);
	std::string DesLocfilename = partfilename + std::string{ "DesLoc" };
	std::string Desinfofilename = partfilename + std::string{ "Desinfo" };
	std::string newinterconfilename = partfilename + std::string{ "newintercon" };
	std::string newexterconfilename = partfilename + std::string{ "newextercon" };
	std::string intergeofilename = partfilename + std::string{ "intergeo" };
	std::string extergeofilename = partfilename + std::string{ "extergeo" };

	{
		std::ofstream file;
		file.open(DesLocfilename, std::ios::binary);
		if (!file.is_open()) {
			std::cerr << "无法打开文件：" << DesLocfilename << std::endl;
			return;
		}

		// 写入向量大小
		size_t size = DesLoc.size();
		file.write(reinterpret_cast<const char*>(&size), sizeof(size));
		// 写入向量数据
		file.write(reinterpret_cast<const char*>(DesLoc.data()), size * sizeof(uint));

		file.close();
	}

	{
		std::ofstream file;
		file.open(Desinfofilename, std::ios::binary);
		if (!file.is_open()) {
			std::cerr << "无法打开文件：" << Desinfofilename << std::endl;
			return;
		}

		// 写入向量大小
		size_t size = Desinfo.size();
		file.write(reinterpret_cast<const char*>(&size), sizeof(size));
		// 写入向量数据
		file.write(reinterpret_cast<const char*>(Desinfo.data()), size * sizeof(uint));

		file.close();
	}
	{
		std::ofstream file;
		file.open(newinterconfilename, std::ios::binary);
		if (!file.is_open()) {
			std::cerr << "无法打开文件：" << newinterconfilename << std::endl;
			return;
		}

		// 写入向量大小
		size_t size = newintercon.size();
		file.write(reinterpret_cast<const char*>(&size), sizeof(size));
		// 写入向量数据
		file.write(reinterpret_cast<const char*>(newintercon.data()), size * sizeof(uint));

		file.close();
	}
	{
		std::ofstream file;
		file.open(newexterconfilename, std::ios::binary);
		if (!file.is_open()) {
			std::cerr << "无法打开文件：" << newexterconfilename << std::endl;
			return;
		}

		// 写入向量大小
		size_t size = newextercon.size();
		file.write(reinterpret_cast<const char*>(&size), sizeof(size));
		// 写入向量数据
		file.write(reinterpret_cast<const char*>(newextercon.data()), size * sizeof(uint));

		file.close();
	}
	{
		std::ofstream file;
		file.open(intergeofilename, std::ios::binary);
		if (!file.is_open()) {
			std::cerr << "无法打开文件：" << intergeofilename << std::endl;
			return;
		}

		// 写入向量大小
		size_t size = intergeo.size();
		file.write(reinterpret_cast<const char*>(&size), sizeof(size));
		// 写入向量数据
		file.write(reinterpret_cast<const char*>(intergeo.data()), size * sizeof(float));

		file.close();
	}
	{
		std::ofstream file;
		file.open(extergeofilename, std::ios::binary);
		if (!file.is_open()) {
			std::cerr << "无法打开文件：" << extergeofilename << std::endl;
			return;
		}

		// 写入向量大小
		size_t size = extergeo.size();
		file.write(reinterpret_cast<const char*>(&size), sizeof(size));
		// 写入向量数据
		file.write(reinterpret_cast<const char*>(extergeo.data()), size * sizeof(float));

		file.close();
	}
}

void Box2d::refresh(float yvalue, float zvalue)
{
	if (yvalue < miny)
		miny = yvalue;
	if (yvalue > maxy)
		maxy = yvalue;

	if (zvalue < minz)
		minz = zvalue;
	if (zvalue > maxz)
		maxz = zvalue;
	return;
}

float Box2d::getSumLength()
{
	return (maxy - miny) + (maxz - minz);
}

void UnitBox::refresh(float xvalue, float yvalue, float zvalue)
{
	if (xvalue < minx)
		minx = xvalue;
	if (xvalue > maxx)
		maxx = xvalue;

	if (yvalue < miny)
		miny = yvalue;
	if (yvalue > maxy)
		maxy = yvalue;

	if (zvalue < minz)
		minz = zvalue;
	if (zvalue > maxz)
		maxz = zvalue;
	return;
}
