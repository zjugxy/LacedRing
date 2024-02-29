#include "NewLWGenerator.h"
#define EMPTYWIRE 0xFFFFFFFF
#define WIREEND 0xFF
#define REVERSEFLAG 0x80000000
#define GREATERSCALEELEM 1.015
#define SCALEXELEMMIN (1e-10)
#define MINREQUIRE 1e-15

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

	//VertexQuantization(mesh);
	NewVertexQuantization(mesh);


	GEObitflowPack();

	PackSimpleLaceWire(mesh);

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

		bool usecompressgeo = true;

		//internal pnts part
			if (usecompressgeo) {
				auto& tempgeo = packinter[mid];
				if (tempgeo.needtransform == false) {
					assert(targets[mid].InternalLW.vertex[0] == EMPTYWIRE);
					//for (int i = 0; i < tempgeo.geobitflow.size(); i += 3) {
					//	float x = *reinterpret_cast<float*>(&tempgeo.geobitflow[i]);
					//	float y = *reinterpret_cast<float*>(&tempgeo.geobitflow[i + 1]);
					//	float z = *reinterpret_cast<float*>(&tempgeo.geobitflow[i + 2]);
					//	geoinfo[vertexbegin + vcnt] = vec4{ x,y,z,1.0 };
					//	vcnt++;
					//}
				}
				else {

					Eigen::MatrixXd dequanmatrix;
					Eigen::Vector3d tranvec;
					ParseTempGeo(tempgeo, dequanmatrix, tranvec);
					uint pntlength = tempgeo.xnum + tempgeo.ynum + tempgeo.znum;

					for (int i = 0; i < targets[mid].InternalLW.vertex.size(); ++i) {
						Eigen::Vector3d  temppnt = ReadData(tempgeo.geobitflow, 12 * 8 + pntlength * i, tempgeo.xnum, tempgeo.ynum, tempgeo.znum);
						auto pnt = dequanmatrix * temppnt + tranvec;
						//auto cmpvalue = mesh.point(mesh.vertex_handle(targets[mid].InternalLW.vertex[i]));
						//Eigen::Vector3d cmpvec{ cmpvalue[0],cmpvalue[1],cmpvalue[2] };
						geoinfo[vertexbegin + vcnt] = vec4{ pnt[0],pnt[1],pnt[2],1.0 };
						vcnt++;
					}
				}
			}
			else {
				for (auto& ver : meshlet.InternalLW.vertex) {
					if (ver == EMPTYWIRE)
						break;
					auto vh = mesh.vertex_handle(ver);
					auto pnt = mesh.point(vh);
					geoinfo[vertexbegin + vcnt] = vec4{ pnt[0],pnt[1],pnt[2],1.0 };
					vcnt++;
				}
			}


		//external pnts part
			//这里的eid,idx不一样写的很大辩
			if (usecompressgeo) {
				for (uint idx = 0; idx < meshlet.ewireidx.size(); ++idx) {
					uint eid = meshlet.ewireidx[idx];
					auto& tempgeo = packexter[eid];
					assert(tempgeo.needtransform == true);

					Eigen::MatrixXd dequanmatrix;
					Eigen::Vector3d tranvec;
					ParseTempGeo(tempgeo, dequanmatrix, tranvec);
					uint pntlength = tempgeo.xnum + tempgeo.ynum + tempgeo.znum;

					if (meshlet.reverse[idx] == false) {
						for (int i = 0; i <gloEwires[eid].vertex.size() - 1; ++i) {
							Eigen::Vector3d  temppnt = ReadData(tempgeo.geobitflow, 12 * 8 + pntlength * i, tempgeo.xnum, tempgeo.ynum, tempgeo.znum);
							auto pnt = dequanmatrix * temppnt + tranvec;
							//auto cmpvalue = mesh.point(mesh.vertex_handle(targets[mid].InternalLW.vertex[i]));
							//Eigen::Vector3d cmpvec{ cmpvalue[0],cmpvalue[1],cmpvalue[2] };
							geoinfo[vertexbegin + vcnt] = vec4{ pnt[0],pnt[1],pnt[2],1.0 };
							vcnt++;
						}
					}
					else {
						for (int i = gloEwires[eid].vertex.size() - 1; i != 0; --i) {
							Eigen::Vector3d  temppnt = ReadData(tempgeo.geobitflow, 12 * 8 + pntlength * i, tempgeo.xnum, tempgeo.ynum, tempgeo.znum);
							auto pnt = dequanmatrix * temppnt + tranvec;
							//auto cmpvalue = mesh.point(mesh.vertex_handle(targets[mid].InternalLW.vertex[i]));
							//Eigen::Vector3d cmpvec{ cmpvalue[0],cmpvalue[1],cmpvalue[2] };
							geoinfo[vertexbegin + vcnt] = vec4{ pnt[0],pnt[1],pnt[2],1.0 };
							vcnt++;
						}
					}

				}

			}
			else {
				for (uint idx = 0; idx < meshlet.ewireidx.size(); ++idx) {
					auto& ewire = gloEwires[meshlet.ewireidx[idx]];
					if (meshlet.reverse[idx] == false) {
						for (uint i = 0; i < ewire.vertex.size() - 1; ++i) {
							auto vh = mesh.vertex_handle(ewire.vertex[i]);
							auto pnt = mesh.point(vh);
							geoinfo[vertexbegin + vcnt] = vec4{ pnt[0],pnt[1],pnt[2],1.0 };
							vcnt++;
						}
					}
					else {
						for (uint i = ewire.vertex.size() - 1; i != 0; --i) {
							auto vh = mesh.vertex_handle(ewire.vertex[i]);
							auto pnt = mesh.point(vh);
							geoinfo[vertexbegin + vcnt] = vec4{ pnt[0],pnt[1],pnt[2],1.0 };
							vcnt++;
						}
					}
				}
			}
			
			assert(vcnt == vidxmap.size());
			
		//}





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
			//PackNoPlane(packexter[eid], mesh, gloEwires[eid]);
			//continue;

			assert(gloEwires[eid].vertex.size() == 2);
			PackExtraPnts(packexter[eid], mesh, gloEwires[eid].vertex);

		}

		auto& tempgeo = packexter[eid];
		auto& ewire = gloEwires[eid];


		Eigen::Vector3d centroid, normal,newcentroid;
		Eigen::Vector3d finalyaxis, finalzaxis;
		Box2d globox;
		UnitBox scalebox;

		if (tempgeo.NeedExtraPnt == false) {
			fitPlane(mesh, ewire.vertex, centroid, normal);
			globox = getYZAxis(mesh, ewire.vertex, normal, finalyaxis, finalzaxis, 15.0 / 180.0 * PI, centroid);
			newcentroid = centroid + (globox.maxy + globox.miny) / 2.0f * finalyaxis + (globox.maxz + globox.minz) / 2.0f * finalzaxis;
			getScaleBox(mesh, ewire.vertex, newcentroid, scalebox, normal, finalyaxis, finalzaxis);
		}
		else {
			//需要额外的点构成平面
			fitPlane(mesh, tempgeo.extraVertex, centroid, normal);
			globox = getYZAxis(mesh, tempgeo.extraVertex, normal, finalyaxis, finalzaxis, 15.0 / 180.0 * PI, centroid);
			newcentroid = centroid + (globox.maxy + globox.miny) / 2.0f * finalyaxis + (globox.maxz + globox.minz) / 2.0f * finalzaxis;
			getScaleBox(mesh, tempgeo.extraVertex, newcentroid, scalebox, normal, finalyaxis, finalzaxis);
		}
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
		//NewSimpleCheck(tempgeo, mesh, ewire.vertex);

		if (firstadd) {
			MeshTransBox = UnitBox{ float(newcentroid[0]),float(newcentroid[1]),float(newcentroid[2])
				,float(newcentroid[0]),float(newcentroid[1]),float(newcentroid[2]) };
			firstadd = false;
		}
		else 
			MeshTransBox.refresh(newcentroid[0], newcentroid[1], newcentroid[2]);
		
	}

	for (uint iid = 0; iid < packinter.size(); ++iid) {
		//个数小于3的就不处理
		if (targets[iid].InternalLW.vertex.size() < 3) {
			//PackNoPlane(packinter[iid], mesh, targets[iid].InternalLW);
			//continue;
			if (targets[iid].InternalLW.vertex[0] == EMPTYWIRE) {
				packinter[iid].needtransform = false;
				continue;
			}
			PackExtraPnts(packinter[iid], mesh, targets[iid].InternalLW.vertex);
		
		}

		auto& tempgeo = packinter[iid];
		auto& ewire = targets[iid].InternalLW;

		Eigen::Vector3d centroid, normal, newcentroid;
		Eigen::Vector3d finalyaxis, finalzaxis;
		Box2d globox;
		UnitBox scalebox;
		if (tempgeo.NeedExtraPnt == false) {
			fitPlane(mesh, ewire.vertex, centroid, normal);
			globox = getYZAxis(mesh, ewire.vertex, normal, finalyaxis, finalzaxis, 15.0 / 180.0 * PI, centroid);
			newcentroid = centroid + (globox.maxy + globox.miny) / 2.0f * finalyaxis + (globox.maxz + globox.minz) / 2.0f * finalzaxis;
			getScaleBox(mesh, ewire.vertex, newcentroid, scalebox, normal, finalyaxis, finalzaxis);
		}
		else {
			//需要额外的点构成平面
			fitPlane(mesh, tempgeo.extraVertex, centroid, normal);
			globox = getYZAxis(mesh, tempgeo.extraVertex, normal, finalyaxis, finalzaxis, 15.0 / 180.0 * PI, centroid);
			newcentroid = centroid + (globox.maxy + globox.miny) / 2.0f * finalyaxis + (globox.maxz + globox.minz) / 2.0f * finalzaxis;
			getScaleBox(mesh, tempgeo.extraVertex, newcentroid, scalebox, normal, finalyaxis, finalzaxis);
		}
		//pack tempgeo
		tempgeo.needtransform = true;
		tempgeo.translation = newcentroid;
		tempgeo.scaleelem = Eigen::Vector3d{ scalebox.maxx,scalebox.maxy,scalebox.maxz };
		Eigen::Matrix3d rotationMatrix;
		rotationMatrix << normal, finalyaxis, finalzaxis;
		Eigen::Vector3d eulerXYZ = rotationMatrix.eulerAngles(0, 1, 2);
		tempgeo.euler = eulerXYZ;

		//NewSimpleCheck(tempgeo, mesh, ewire.vertex);

		MeshTransBox.refresh(newcentroid[0], newcentroid[1], newcentroid[2]);
	}
	

	//使用8-bit将 eulerxyz --> 映射到对应的 0-2PI值域
	UcharRadGen();

	//for (uint eid = 0; eid < gloEwires.size(); ++eid)
	//	if(packexter[eid].needtransform == true)
	//		NewSimpleCheck(packexter[eid], mesh, gloEwires[eid].vertex);

	//translation直接确定好一点，就按照原来的进行量化
	TranslatePack(mesh, MeshTransBox);

	FinalScaleGen(mesh);
	//pack data to compress;
	CheckByDequantize(mesh);


	BitSortGen(30,2);
	VertexBitGen(mesh,0.001);

	//pack into bit flow


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

void NewLWGenerator::BitSortGen(int highestnum = 30, int minvalue = 2)
{
	for (int i = minvalue*3; i < highestnum; ++i) {

		for (int a = minvalue; a < i; ++a)
			for (int b = minvalue; b < i - a; ++b) {
				int c = i - a - b;
				if (c < minvalue)continue;
				std::vector<int> temp{ a,b,c };
				bitnums.push_back(temp);
			}

	}

}

void NewLWGenerator::VertexBitGen(const MyMesh& mesh,float errorpercent)
{
	//使用 meshtranslation中的数据来作为一个标准
	float xerror = (MeshTransBox.maxx - MeshTransBox.minx) * errorpercent;
	float yerror = (MeshTransBox.maxy - MeshTransBox.miny) * errorpercent;
	float zerror = (MeshTransBox.maxz - MeshTransBox.minz) * errorpercent;

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
		Eigen::Matrix3d rotationMatrix = (rotationX.matrix() * rotationY.matrix() * rotationZ.matrix());
		Eigen::Vector3d tranvec = {
			MeshTransBox.minx + (tempgeo.translatex / 255.0) * (MeshTransBox.maxx - MeshTransBox.minx),
			MeshTransBox.miny + (tempgeo.translatey / 255.0) * (MeshTransBox.maxy - MeshTransBox.miny),
			MeshTransBox.minz + (tempgeo.translatez / 255.0) * (MeshTransBox.maxz - MeshTransBox.minz)
		};
		float descalex = MeshScaleBox.minx * std::pow(MeshScaleBox.maxx / MeshScaleBox.minx, tempgeo.scalex / 255.0f);
		float descaley = MeshScaleBox.miny * std::pow(MeshScaleBox.maxy / MeshScaleBox.miny, tempgeo.scaley / 255.0f);
		float descalez = MeshScaleBox.minz * std::pow(MeshScaleBox.maxz / MeshScaleBox.minz, tempgeo.scalez / 255.0f);

		Eigen::DiagonalMatrix<double, 3> scalema;
		scalema.diagonal() << descalex, descaley, descalez;
		Eigen::MatrixXd dequanmatrix = rotationMatrix * scalema;

		for (int type3 = 0; type3 < bitnums.size(); ++type3) {
			int xnum = bitnums[type3][0];
			int ynum = bitnums[type3][1];
			int znum = bitnums[type3][2];

			bool flag = true;

			for (uint vid = 0; vid < gloEwires[eid].vertex.size(); ++vid) {
				auto ver = gloEwires[eid].vertex[vid];
				auto pnt = mesh.point(mesh.vertex_handle(ver));
				auto rawdata = Eigen::Vector3d{ pnt[0],pnt[1],pnt[2] };
				Eigen::Vector3d compressedata = tempgeo.datatocompress[vid];

				Eigen::Vector3d tempdata = cutfloatByBit(compressedata, xnum, ynum, znum);
				Eigen::Vector3d cmpdata = dequanmatrix * tempdata + tranvec;
				if ((abs(cmpdata.x() - rawdata.x()) > xerror) || (abs(cmpdata.y() - rawdata.y()) > yerror) || (abs(cmpdata.z() - rawdata.z()) > zerror)) {
					flag = false; break;
				}
			}
			if (flag == true) {
				tempgeo.xnum = xnum; tempgeo.ynum = ynum; tempgeo.znum = znum;
				break;
			}
			assert(type3 <(bitnums.size()-1) && "error in vertex bit gen");
		}

		std::cout << "eid is " <<eid << " " << uint(tempgeo.xnum) << " " << uint(tempgeo.ynum) << " " << uint(tempgeo.znum) << std::endl;

#ifndef NDEBUG
		//test
		{
			int xnum = tempgeo.xnum, ynum = tempgeo.ynum, znum = tempgeo.znum;
			for (uint vid = 0; vid < gloEwires[eid].vertex.size(); ++vid) {
				auto ver = gloEwires[eid].vertex[vid];
				auto pnt = mesh.point(mesh.vertex_handle(ver));
				auto rawdata = Eigen::Vector3d{ pnt[0],pnt[1],pnt[2] };
				Eigen::Vector3d compressedata = tempgeo.datatocompress[vid];

				Eigen::Vector3d tempdata = cutfloatByBit(compressedata, xnum, ynum, znum);
				Eigen::Vector3d cmpdata = dequanmatrix * tempdata + tranvec;
				std::cout << (cmpdata - rawdata).norm() << std::endl;
			}
		}
#endif // !NDEBUG


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
		Eigen::Matrix3d rotationMatrix = (rotationX.matrix() * rotationY.matrix() * rotationZ.matrix());
		Eigen::Vector3d tranvec = {
			MeshTransBox.minx + (tempgeo.translatex / 255.0) * (MeshTransBox.maxx - MeshTransBox.minx),
			MeshTransBox.miny + (tempgeo.translatey / 255.0) * (MeshTransBox.maxy - MeshTransBox.miny),
			MeshTransBox.minz + (tempgeo.translatez / 255.0) * (MeshTransBox.maxz - MeshTransBox.minz)
		};
		float descalex = MeshScaleBox.minx * std::pow(MeshScaleBox.maxx / MeshScaleBox.minx, tempgeo.scalex / 255.0f);
		float descaley = MeshScaleBox.miny * std::pow(MeshScaleBox.maxy / MeshScaleBox.miny, tempgeo.scaley / 255.0f);
		float descalez = MeshScaleBox.minz * std::pow(MeshScaleBox.maxz / MeshScaleBox.minz, tempgeo.scalez / 255.0f);

		Eigen::DiagonalMatrix<double, 3> scalema;
		scalema.diagonal() << descalex, descaley, descalez;
		Eigen::MatrixXd dequanmatrix = rotationMatrix * scalema;

		for (int type3 = 0; type3 < bitnums.size(); ++type3) {
			int xnum = bitnums[type3][0];
			int ynum = bitnums[type3][1];
			int znum = bitnums[type3][2];

			bool flag = true;

			for (uint vid = 0; vid < targets[iid].InternalLW.vertex.size(); ++vid) {
				auto ver = targets[iid].InternalLW.vertex[vid];
				auto pnt = mesh.point(mesh.vertex_handle(ver));
				auto rawdata = Eigen::Vector3d{ pnt[0],pnt[1],pnt[2] };
				Eigen::Vector3d compressedata = tempgeo.datatocompress[vid];

				Eigen::Vector3d tempdata = cutfloatByBit(compressedata, xnum, ynum, znum);
				Eigen::Vector3d cmpdata = dequanmatrix * tempdata + tranvec;
				if ((abs(cmpdata.x() - rawdata.x()) > xerror) || (abs(cmpdata.y() - rawdata.y()) > yerror) || (abs(cmpdata.z() - rawdata.z()) > zerror)) {
					flag = false; break;
				}
			}
			if (flag == true) {
				tempgeo.xnum = xnum; tempgeo.ynum = ynum; tempgeo.znum = znum;
				break;
			}
			assert(type3 < (bitnums.size() - 1) && "error in vertex bit gen");

		}

		std::cout << "iid is " << iid << " " << uint(tempgeo.xnum) << " " << uint(tempgeo.ynum) << " " << uint(tempgeo.znum) << std::endl;
#ifndef NDEBUG
		//test
		{
			int xnum = tempgeo.xnum, ynum = tempgeo.ynum, znum = tempgeo.znum;
			for (uint vid = 0; vid < targets[iid].InternalLW.vertex.size(); ++vid) {
				auto ver = targets[iid].InternalLW.vertex[vid];
				auto pnt = mesh.point(mesh.vertex_handle(ver));
				auto rawdata = Eigen::Vector3d{ pnt[0],pnt[1],pnt[2] };
				Eigen::Vector3d compressedata = tempgeo.datatocompress[vid];

				Eigen::Vector3d tempdata = cutfloatByBit(compressedata, xnum, ynum, znum);
				Eigen::Vector3d cmpdata = dequanmatrix * tempdata + tranvec;
				std::cout << (cmpdata - rawdata).norm() << std::endl;
			}
		}
#endif // !NDEBUG
	}

	return;
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

void NewLWGenerator::NewSimpleCheck(PackGEO& tempgeo, const MyMesh& mesh, const std::vector<uint>& verset)
{
	Eigen::Vector3d eulerXYZ = tempgeo.euler;
	Eigen::AngleAxisd rotationX(eulerXYZ[0], Eigen::Vector3d::UnitX());
	Eigen::AngleAxisd rotationY(eulerXYZ[1], Eigen::Vector3d::UnitY());
	Eigen::AngleAxisd rotationZ(eulerXYZ[2], Eigen::Vector3d::UnitZ());
	//dequantize
	Eigen::Matrix3d rotationMatrix = rotationX.matrix() * rotationY.matrix() * rotationZ.matrix();
	Eigen::DiagonalMatrix<double, 3> scaleMatirx;
	if(tempgeo.isXPlatForm == false)
		scaleMatirx.diagonal() << tempgeo.scaleelem[0], tempgeo.scaleelem[1], tempgeo.scaleelem[2];
	else
		scaleMatirx.diagonal() << 0.0, tempgeo.scaleelem[1], tempgeo.scaleelem[2];

	//quanitize
	Eigen::Matrix3d revrotat = rotationMatrix.transpose();
	Eigen::DiagonalMatrix<double, 3> revscale;
	if(tempgeo.isXPlatForm == false)
		revscale.diagonal() << 1.0/tempgeo.scaleelem[0], 1.0/tempgeo.scaleelem[1], 1.0/tempgeo.scaleelem[2];
	else
		revscale.diagonal() << 0.0, 1.0 / tempgeo.scaleelem[1], 1.0 / tempgeo.scaleelem[2];


	Eigen::Matrix3d com = revscale * revrotat ;
	Eigen::Matrix3d decom =   rotationMatrix* scaleMatirx;
	//
	float xrecord = 0, yrecord = 0, zrecord = 0;

	for (const auto& ver : verset) {
		auto vh = mesh.vertex_handle(ver);
		auto pnt = mesh.point(vh);
		auto rawdata = Eigen::Vector3d{ pnt[0],pnt[1],pnt[2] };

		Eigen::Vector3d unitdata = com * (rawdata - tempgeo.translation);

		if (abs(unitdata.x() > xrecord))xrecord = abs(unitdata.x());
		if (abs(unitdata.y() > yrecord))yrecord = abs(unitdata.y());
		if (abs(unitdata.z() > zrecord))zrecord = abs(unitdata.z());


		tempgeo.dataunit.push_back(unitdata);
		Eigen::Vector3d cmpdata = decom * unitdata + tempgeo.translation;

		double diff = (cmpdata - rawdata).norm();
		assert(diff < 1e-10);
		//std::cout << "distance cmp and raw data" << diff << std::endl;
	}

	if ((xrecord > 1.00001) || (yrecord > 1.00001) || (zrecord > 1.00001))
		std::cout << (tempgeo.scaleelem.x()<1e-10) << " " << xrecord<<" "<<yrecord << " " << zrecord << std::endl;

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
		if (tempgeo.isXPlatForm == true)descalex = 0.0;


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
			auto debugvalue = result;
			LimitEigen(result);
			tempgeo.datatocompress.push_back(result);
			Eigen::Vector3d dequandata = dequanmatrix * result + tranvec;

			//assert(result.x() < 1.001);
			//assert(result.y() < 1.001);
			//assert(result.z() < 1.001);

			//由于scale每一个元素通过8bit表示 1符号位7数值位
			std::cout<<"差距是 " << (dequandata - rawdata).norm() << std::endl;
			if ((dequandata - rawdata).norm() > MINREQUIRE) {
				std::cout << "debug 需要放大scale elem的放大率" << std::endl;
				assert(false);
			}
		}


	}

	for (uint iid = 0; iid < packinter.size(); ++iid) {
		auto& tempgeo = packinter[iid];
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
		if (tempgeo.isXPlatForm == true)descalex = 0.0;


		Eigen::DiagonalMatrix<double, 3> revscale, scalema;
		revscale.diagonal() << 1.0 / descalex, 1.0 / descaley, 1.0 / descalez;
		scalema.diagonal() << descalex, descaley, descalez;

		Eigen::MatrixXd compressmatrix = revscale * rotationMatrix;
		Eigen::MatrixXd dequanmatrix = rotationMatrix.transpose() * scalema;

		auto vertexset = targets[iid].InternalLW.vertex;
		for (auto ver : vertexset) {
			auto pnt = mesh.point(mesh.vertex_handle(ver));
			auto rawdata = Eigen::Vector3d{ pnt[0],pnt[1],pnt[2] };

			Eigen::Vector3d result = compressmatrix * (rawdata - tranvec);
			auto debugvalue = result;
			LimitEigen(result);
			tempgeo.datatocompress.push_back(result);
			Eigen::Vector3d dequandata = dequanmatrix * result + tranvec;

			//assert(result.x() < 1.001);
			//assert(result.y() < 1.001);
			//assert(result.z() < 1.001);

			//由于scale每一个元素通过8bit表示 1符号位7数值位
			std::cout << "差距是 " << (dequandata - rawdata).norm() << std::endl;
			if ((dequandata - rawdata).norm() > MINREQUIRE) {
				std::cout << "debug 需要放大scale elem的放大率" << std::endl;
				assert(false);
			}
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
		if (tempgeo.NeedExtraPnt)
			vertexset = tempgeo.extraVertex;

		float tempx = 0, tempy = 0,tempz = 0;
		for (auto ver : vertexset) {
			auto pnt = mesh.point(mesh.vertex_handle(ver));
			Eigen::Vector3d point{ pnt[0],pnt[1],pnt[2] };
			point -= tranvec;
			Eigen::Vector3d scalevec = rotationMatrix * point;
			if (abs(scalevec.x()) > tempx)tempx = abs(scalevec.x());
			if (abs(scalevec.y()) > tempy)tempy = abs(scalevec.y());
			if (abs(scalevec.z()) > tempz)tempz = abs(scalevec.z());
		}
		xvalues.push_back(tempx*GREATERSCALEELEM); yvalues.push_back(tempy * GREATERSCALEELEM); zvalues.push_back(tempz * GREATERSCALEELEM);
		//std::cout << "ewire " << eid << " " << tempx << " " << tempy << " " << tempz << std::endl;

		//如何 scale elem过于小的话使用指数表示有一定精度问题
		assert((tempx > SCALEXELEMMIN) && (tempy > SCALEXELEMMIN) && (tempz > SCALEXELEMMIN));
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
		if (tempgeo.NeedExtraPnt)
			vertexset = tempgeo.extraVertex;

		float tempx = 0, tempy = 0, tempz = 0;
		for (auto ver : vertexset) {
			auto pnt = mesh.point(mesh.vertex_handle(ver));
			Eigen::Vector3d point{ pnt[0],pnt[1],pnt[2] };
			point -= tranvec;
			Eigen::Vector3d scalevec = rotationMatrix * point;
			if (abs(scalevec.x()) > tempx)tempx = abs(scalevec.x());
			if (abs(scalevec.y()) > tempy)tempy = abs(scalevec.y());
			if (abs(scalevec.z()) > tempz)tempz = abs(scalevec.z());
		}
		xvalues.push_back(tempx * GREATERSCALEELEM); yvalues.push_back(tempy * GREATERSCALEELEM); zvalues.push_back(tempz * GREATERSCALEELEM);
		//std::cout << "iwire " << iid << " " << tempx << " " << tempy << " " << tempz << std::endl;
		assert((tempx > SCALEXELEMMIN) && (tempy > SCALEXELEMMIN) && (tempz > SCALEXELEMMIN));

	}


	//由于引入了误差的原因，反而问题没那么大了
	MeshScaleBox = UnitBox{ xvalues.front(),yvalues.front(),zvalues.front(),xvalues.front(),yvalues.front(),zvalues.front() };
	//pack 部分
	for (auto elem : xvalues) {
		if (elem > MeshScaleBox.maxx)MeshScaleBox.maxx = elem;
	}
	for (auto elem : yvalues) {
		if (elem > MeshScaleBox.maxy)MeshScaleBox.maxy = elem;
	}
	for (auto elem : zvalues) {
		if (elem > MeshScaleBox.maxz)MeshScaleBox.maxz = elem;
	}

	uint meshletid = 0;

	float undernumx = std::log(MeshScaleBox.maxx / MeshScaleBox.minx);
	float undernumy = std::log(MeshScaleBox.maxy / MeshScaleBox.miny);
	float undernumz = std::log(MeshScaleBox.maxz / MeshScaleBox.minz);


	for (uint eid = 0; eid < packexter.size(); ++eid) {
		auto& tempgeo = packexter[eid];
		if (tempgeo.needtransform == false)continue;
		tempgeo.scalex = myclamp(255 * (std::log(xvalues[meshletid] / MeshScaleBox.minx) / undernumx),0,255);
		tempgeo.scaley = myclamp(255 * (std::log(yvalues[meshletid] / MeshScaleBox.miny) / undernumy),0,255);
		tempgeo.scalez = myclamp(255 * (std::log(zvalues[meshletid] / MeshScaleBox.minz) / undernumz),0,255);
		meshletid++;
	}

	for (uint iid = 0; iid < packinter.size(); ++iid) {
		auto& tempgeo = packinter[iid];
		if (tempgeo.needtransform == false)continue;
		tempgeo.scalex = myclamp(255 * (std::log(xvalues[meshletid] / MeshScaleBox.minx) / undernumx), 0, 255);
		tempgeo.scaley = myclamp(255 * (std::log(yvalues[meshletid] / MeshScaleBox.miny) / undernumy), 0, 255);
		tempgeo.scalez = myclamp(255 * (std::log(zvalues[meshletid] / MeshScaleBox.minz) / undernumz), 0, 255);
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
	int cmpvalue = std::round(value);

	if (cmpvalue > highbound)return highbound;
	else if (cmpvalue < lowbound)return lowbound;
	return uchar(cmpvalue);
}

void NewLWGenerator::LimitEigen(Eigen::Vector3d& vec)
{
	if (vec.x() < -1.0)vec.x() = -1.0;
	if (vec.x() >  1.0)vec.x() =  1.0;

	if (vec.y() < -1.0)vec.y() = -1.0;
	if (vec.y() > 1.0)vec.y() = 1.0;

	if (vec.z() < -1.0)vec.z() = -1.0;
	if (vec.z() > 1.0)vec.z() = 1.0;
}

// float*3 在该精度之下具体会变成什么值
// -1.0 --> 1.0  map to  0 --> 2^num-1
// 0.0 --> 2.0 map to  0 --> 2^num-1



Eigen::Vector3d NewLWGenerator::cutfloatByBit(Eigen::Vector3d rawvalue, uint xnum, uint ynum, uint znum)
{
	uint xhigh = (1 << xnum) - 1;
	uint yhigh = (1 << ynum) - 1;
	uint zhigh = (1 << znum) - 1;

	float x = std::round((rawvalue.x()+1)/2.0 * xhigh);//value 对应uint值
	float y = std::round((rawvalue.y()+1)/2.0 * yhigh);
	float z = std::round((rawvalue.z()+1)/2.0 * zhigh);
	auto result = Eigen::Vector3d{
		x/xhigh*2.0 -1.0,y/yhigh*2.0 -1.0,z/zhigh*2.0-1.0
	};
	return result;
}

void NewLWGenerator::GEObitflowPack()
{
	for (uint eid = 0; eid < packexter.size(); ++eid) {
		auto& tempgeo = packexter[eid];
		
		//<3
		if (tempgeo.needtransform == false) {
			//代码修改之后这部分 ewire needtranform总是 true
			assert(false);
			//assert(gloEwires[eid].vertex.size() < 3);

			//for (auto pnt : tempgeo.nopackgeo) {
			//	uint uintvalue = *reinterpret_cast<uint*>(&pnt.x);
			//	tempgeo.geobitflow.push_back(uintvalue);
			//	uintvalue = *reinterpret_cast<uint*>(&pnt.y);
			//	tempgeo.geobitflow.push_back(uintvalue);
			//	uintvalue = *reinterpret_cast<uint*>(&pnt.z);
			//	tempgeo.geobitflow.push_back(uintvalue);
			//}
		}
		else {
			tempgeo.geobitflow.push_back(PackChar4Uint(tempgeo.rotatx, tempgeo.rotaty, tempgeo.rotatz, tempgeo.translatex));
			tempgeo.geobitflow.push_back(PackChar4Uint(tempgeo.translatey, tempgeo.translatez, tempgeo.scalex, tempgeo.scaley));
			tempgeo.geobitflow.push_back(PackChar4Uint(tempgeo.scalez, tempgeo.xnum, tempgeo.ynum, tempgeo.znum));

			uint xnum = tempgeo.xnum,ynum = tempgeo.ynum,znum = tempgeo.znum;
			uint pntlength = xnum + ynum + znum;
			tempgeo.geobitflow.resize(3 + (31 + pntlength * tempgeo.datatocompress.size()) / 32);
			uint startbitidx = 3 * 32;

			for (auto pnt : tempgeo.datatocompress) {
				InsertGeoValue(tempgeo.geobitflow, startbitidx, xnum, pnt.x());
				InsertGeoValue(tempgeo.geobitflow, startbitidx+xnum, ynum, pnt.y());
				InsertGeoValue(tempgeo.geobitflow, startbitidx+xnum+ynum, znum, pnt.z());
				startbitidx += pntlength;
			}
		}
	
	
	
	}

	for (uint iid = 0; iid < packinter.size(); ++iid) {
		auto& tempgeo = packinter[iid];

		//<3
		if (tempgeo.needtransform == false) {
			assert(targets[iid].InternalLW.vertex[0]==EMPTYWIRE);

			//for (auto pnt : tempgeo.nopackgeo) {
			//	uint uintvalue = *reinterpret_cast<uint*>(&pnt.x);
			//	tempgeo.geobitflow.push_back(uintvalue);
			//	uintvalue = *reinterpret_cast<uint*>(&pnt.y);
			//	tempgeo.geobitflow.push_back(uintvalue);
			//	uintvalue = *reinterpret_cast<uint*>(&pnt.z);
			//	tempgeo.geobitflow.push_back(uintvalue);
			//}
		}
		else {
			tempgeo.geobitflow.push_back(PackChar4Uint(tempgeo.rotatx, tempgeo.rotaty, tempgeo.rotatz, tempgeo.translatex));
			tempgeo.geobitflow.push_back(PackChar4Uint(tempgeo.translatey, tempgeo.translatez, tempgeo.scalex, tempgeo.scaley));
			tempgeo.geobitflow.push_back(PackChar4Uint(tempgeo.scalez, tempgeo.xnum, tempgeo.ynum, tempgeo.znum));

			uint xnum = tempgeo.xnum, ynum = tempgeo.ynum, znum = tempgeo.znum;
			uint pntlength = xnum + ynum + znum;
			tempgeo.geobitflow.resize(3 + (31 + pntlength * tempgeo.datatocompress.size()) / 32);
			uint startbitidx = 3 * 32;

			for (auto pnt : tempgeo.datatocompress) {
				InsertGeoValue(tempgeo.geobitflow, startbitidx, xnum, pnt.x());
				InsertGeoValue(tempgeo.geobitflow, startbitidx + xnum, ynum, pnt.y());
				InsertGeoValue(tempgeo.geobitflow, startbitidx + xnum + ynum, znum, pnt.z());
				startbitidx += pntlength;
			}
		}



	}

}

void NewLWGenerator::InsertGeoValue(std::vector<uint>& bitflow, uint startidx, uint length, float rawvalue)
{
	uint highvalue = (1 << length) - 1;
	assert((rawvalue >= -1.0) && (rawvalue <= 1.0));
	uint cutvalue = std::round((rawvalue+1.0)/2.0 * highvalue);
	// startidx --> startidx+length
	uint endidx = startidx + length;

	//[startidx,endidx)
	//从低位到高位开始覆写
	if ((endidx-1) / 32 == startidx / 32) {
		uint& value = bitflow[startidx / 32];
		uint offset = startidx % 32;
		//printBits(value);
		//printBits(cutvalue);
		value = value | (cutvalue << offset);
		//printBits(value);
	}
	else {
		uint offset = startidx % 32;
		uint firstpart = cutvalue << offset;
		uint secondpart = cutvalue >> (32 - offset);
		//std::cout << "special insert test" << std::endl;
		//printBits(firstpart);
		//printBits(secondpart);
		//printBits(bitflow[startidx / 32]);
		bitflow[startidx / 32] = bitflow[startidx / 32] | firstpart;
		bitflow[startidx / 32 + 1] = bitflow[startidx / 32 + 1] | secondpart;

		//printBits(bitflow[startidx / 32]);
		//printBits(bitflow[startidx / 32+1]);

	}

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

void NewLWGenerator::printBits(uint32_t value) {
	uint32_t mask = 1 << 31;  // 从最高位开始

	for (int i = 0; i < 32; i++) {
		bool bit = (value & mask) != 0;
		std::cout << (bit ? "1" : "0");
		mask >>= 1;
	}

	std::cout << std::endl;
}

Eigen::Vector3d NewLWGenerator::ReadData(const std::vector<uint>& geobits, uint startidx, uint xnum, uint ynum, uint znum)
{
	float x = ReadFloat(geobits, startidx, xnum);
	float y = ReadFloat(geobits, startidx+xnum, ynum);
	float z = ReadFloat(geobits, startidx+xnum+ynum, znum);

	return Eigen::Vector3d{x,y,z};
}

float NewLWGenerator::ReadFloat(const std::vector<uint>& geobits, uint startidx, uint num)
{
	int start = startidx / 32;
	int end = (startidx + num -1) / 32;

	if (start == end) {
		uint value = geobits[start];
		uint offset = startidx % 32;
		uint mask = (1 << (num )) - 1;
		value = value >> offset;
		value = value & mask;

		uint high = (1 << num) - 1;
		float ret = float(value) / high * 2.0 - 1.0;
		return ret;
	}
	else {
		uint offset = startidx % 32;
		uint highmasknum = 32 - offset;
		uint lowmasknum = num - highmasknum;
		uint part1 = geobits[start] >> offset;

		uint lowmask = (1 << (lowmasknum )) - 1;
		uint part2 = geobits[start + 1] & lowmask;
		uint value = part1 + (part2 << highmasknum);

		uint high = (1 << num) - 1;
		float ret = float(value) / high * 2.0 - 1.0;

		if (abs(ret) > 2.0)
			std::cout << "debug" << std::endl;
		return ret;
	}

}

void NewLWGenerator::PackExtraPnts(PackGEO& tempgeo, const MyMesh& mesh, const std::vector<uint>& vertices)
{
	if (vertices[0] == EMPTYWIRE)return;
	//找共享点

	tempgeo.NeedExtraPnt = true;
	if (vertices.size() == 2) {
		uint pnt0 = vertices[0];
		uint pnt1 = vertices[1];

		tempgeo.extraVertex.push_back(pnt0);
		tempgeo.extraVertex.push_back(pnt1);


		std::unordered_set<uint> pntsaround;
		auto vh = mesh.vertex_handle(pnt0);
		for (auto cvvit = mesh.cvv_iter(vh); cvvit.is_valid(); ++cvvit) {
			pntsaround.insert(cvvit->idx());
		}

		vh = mesh.vertex_handle(pnt1);
		for (auto cvvit = mesh.cvv_iter(vh); cvvit.is_valid(); ++cvvit) {
			if (pntsaround.find(cvvit->idx()) != pntsaround.end())
			{
				tempgeo.extraVertex.push_back(cvvit->idx());
				return;
			}
		}
		//internalwire会出现两个不相邻的点
		tempgeo.extraVertex.push_back(*pntsaround.begin());
		assert(tempgeo.extraVertex.size() == 3);
		return;

	}

	if (vertices.size() == 1) {
		uint pnt0 = vertices[0];

		tempgeo.extraVertex.push_back(pnt0);
		auto vh = mesh.vertex_handle(pnt0);
		for (auto cvvit = mesh.cvv_iter(vh); cvvit.is_valid(); ++cvvit) {
			tempgeo.extraVertex.push_back(cvvit->idx());
			if (tempgeo.extraVertex.size() == 3)
				return;
		}
	}
	assert(false);
}

void NewLWGenerator::ParseTempGeo(const PackGEO& tempgeo, Eigen::MatrixXd& dequanmatrix, Eigen::Vector3d& tranvec)
{

	Eigen::Vector3d eulerXYZ{
						float(tempgeo.rotatx) / 255.0 * TWOPI,
						float(tempgeo.rotaty) / 255.0 * TWOPI,
						float(tempgeo.rotatz) / 255.0 * TWOPI
	};
	Eigen::AngleAxisd rotationX(eulerXYZ[0], Eigen::Vector3d::UnitX());
	Eigen::AngleAxisd rotationY(eulerXYZ[1], Eigen::Vector3d::UnitY());
	Eigen::AngleAxisd rotationZ(eulerXYZ[2], Eigen::Vector3d::UnitZ());
	//quantize
	Eigen::Matrix3d rotationMatrix = (rotationX.matrix() * rotationY.matrix() * rotationZ.matrix());
	tranvec = {
		MeshTransBox.minx + (tempgeo.translatex / 255.0) * (MeshTransBox.maxx - MeshTransBox.minx),
		MeshTransBox.miny + (tempgeo.translatey / 255.0) * (MeshTransBox.maxy - MeshTransBox.miny),
		MeshTransBox.minz + (tempgeo.translatez / 255.0) * (MeshTransBox.maxz - MeshTransBox.minz)
	};
	float descalex = MeshScaleBox.minx * std::pow(MeshScaleBox.maxx / MeshScaleBox.minx, tempgeo.scalex / 255.0f);
	float descaley = MeshScaleBox.miny * std::pow(MeshScaleBox.maxy / MeshScaleBox.miny, tempgeo.scaley / 255.0f);
	float descalez = MeshScaleBox.minz * std::pow(MeshScaleBox.maxz / MeshScaleBox.minz, tempgeo.scalez / 255.0f);

	Eigen::DiagonalMatrix<double, 3> scalema;
	scalema.diagonal() << descalex, descaley, descalez;
	dequanmatrix = rotationMatrix * scalema;
}
