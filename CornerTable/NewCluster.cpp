#include "NewCluster.h"
#include <Eigen/Dense>
#include<chrono>

bool CompareScores(const std::pair<uint32_t, float>& a, const std::pair<uint32_t, float>& b);

NewCluster::NewCluster(uint maxverts, uint maxtris, const MyMesh& mesh)
{
	auto startTime = std::chrono::steady_clock::now();

	nfaces = mesh.n_faces();
	nvertices = mesh.n_vertices();

	mymeshlets.clear();
	mymeshlets.emplace_back();
	Meshlet_built* curr = &mymeshlets.back();

	std::vector<std::pair<uint, float>> candidates;
	std::unordered_set<uint> candidateCheck;
	std::vector<bool> checklist(nfaces,false);
	
	maxf = maxtris;
	maxv = maxverts;
	//points
	//normals
	std::vector<vec3> positions;
	std::vector<vec3> normals;


	uint triIndex = 0;
	candidates.push_back(std::make_pair(triIndex, 0.0f));
	candidateCheck.insert(triIndex);

	while (!candidates.empty())
	{
		uint faceindex = candidates.back().first;
		candidates.pop_back();

		auto fh = mesh.face_handle(faceindex);
		std::array<uint, 3> pnts;
		{
			int temp = 0;
			MyMesh::ConstFaceVertexCCWIter fv_it = mesh.cfv_ccwiter(fh);
			for (; fv_it.is_valid(); ++fv_it)
			{
				MyMesh::VertexHandle vertex = *fv_it;
				pnts[temp++] = vertex.idx();
			}
		}

		assert(pnts[0] < nvertices);
		assert(pnts[1] < nvertices);
		assert(pnts[2] < nvertices);
		//有问题啊
		if (AddToMeshlet(maxverts, maxtris, curr, pnts,faceindex)) {
			//faceindex对应三角形被插入meshlet
			checklist[faceindex] = true;

			//计算 center radius normal
			for (uint i = 0; i < 3; ++i) {
				uint vidx = pnts[i];
				auto vh = mesh.vertex_handle(vidx);
				auto geom_vertex = mesh.point(vh);
				positions.emplace_back(geom_vertex[0], geom_vertex[1], geom_vertex[2]);
			}
			
			{
				const vec3& vertex1 = positions[positions.size() - 3];
				const vec3& vertex2 = positions[positions.size() - 2];
				const vec3& vertex3 = positions[positions.size() - 1];

				auto edge1 = vertex2 - vertex1;
				auto edge2 = vertex3 - vertex2;

				auto normal = glm::normalize(glm::cross(edge1, edge2));
				normals.push_back(normal);
			}
			float radius;
			vec3 center;
			CalculateCenterRadius(positions, radius, center);

			float _useless;
			vec3 normal;
			CalculateCenterRadius(normals, _useless, normal);

			//遍历face相邻的三个face 并且更新 candidates
			{
				MyMesh::ConstFaceFaceCCWIter ff_it = mesh.cff_ccwiter(fh);
				for (; ff_it.is_valid(); ++ff_it) {
					auto fidx = ff_it->idx();
					if (fidx<0 || fidx>nfaces)
						continue;
					if (checklist[fidx] == true)
						continue;
					if (candidateCheck.count(fidx))
						continue;

					candidates.push_back(std::make_pair(fidx, FLT_MAX));
					candidateCheck.insert(fidx);
				}
			}

			for (auto& elem : candidates) {
				uint fidx = elem.first;
				std::array<vec3, 3> tri;
				std::array<uint, 3> triidxs;
				auto nfh = mesh.face_handle(fidx);

				int temp = 0;
				MyMesh::ConstFaceVertexCCWIter fv_it = mesh.cfv_ccwiter(nfh);
				for (; fv_it.is_valid(); ++fv_it)
				{
					MyMesh::VertexHandle vertex = *fv_it;
					auto pnt = mesh.point(vertex);
					tri[temp] = vec3{ pnt[0],pnt[1],pnt[2] };
					triidxs[temp++] = vertex.idx();
				}

				elem.second = ComputeScore(curr, center, radius, normal, tri, triidxs);
			}

			if (IsMeshletFull(curr, maxverts, maxtris)||candidates.empty()) {
				positions.clear();
				normals.clear();
				candidateCheck.clear();
				if (!candidates.empty()) {
					candidates[0] = candidates.back();
					candidates.resize(1);
					candidateCheck.insert(candidates[0].first);
				}
				mymeshlets.emplace_back();
				curr = &mymeshlets.back();
			}
			else
			{
				std::sort(candidates.begin(), candidates.end(), CompareScores);
			}
		}
		else {
			if (candidates.empty())
			{
				positions.clear();
				normals.clear();
				candidateCheck.clear();

				mymeshlets.emplace_back();
				curr = &mymeshlets.back();
			}
		}

		if (candidates.empty())
		{
			while (triIndex < nfaces && checklist[triIndex])
				++triIndex;

			if (triIndex == nfaces)
				break;

			candidates.push_back(std::make_pair(triIndex, 0.0f));
			candidateCheck.insert(triIndex);
		}

	}

	if (mymeshlets.back().faces.empty())
	{
		mymeshlets.pop_back();
	}

	//可以在这里做一个模拟退火？
	// cut sharp triangles
	auto endTime = std::chrono::steady_clock::now();
	std::cout << "first cluster down" << std::endl;
	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);
	std::cout << "代码执行时间: " << duration.count() << " 毫秒" << std::endl;

	 startTime = std::chrono::steady_clock::now();
	//识别环条状meshlet以及近似还条状的meshlet,记录并分割
	while (SharpTriCut(mesh)) {};
	std::cout << "sharp cut down" << std::endl;
	 endTime = std::chrono::steady_clock::now();
	std::cout << "first cluster down" << std::endl;
	 duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);
	std::cout << "代码执行时间: " << duration.count() << " 毫秒" << std::endl;

	 startTime = std::chrono::steady_clock::now();
	TryShapeHeal(mesh);
	std::cout << "shape heal down" << std::endl;
	 endTime = std::chrono::steady_clock::now();
	std::cout << "first cluster down" << std::endl;
	 duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);
	std::cout << "代码执行时间: " << duration.count() << " 毫秒" << std::endl;

	PackOldmeshlets();
	GenerateEwire(mesh);

}

bool NewCluster::AddToMeshlet(int maxv, int maxt, Meshlet_built* curr, std::array<uint, 3> pnts,uint& faceid)
{
	Meshlet_built& meshlet = *curr;
	if (meshlet.faces.size() == maxt)
		return false;
	if (meshlet.vertices.size() == maxv)
		return false;


	uint newvertex = 0;
	for (uint i = 0; i < 3u; ++i) {
		if (meshlet.vertices.find(pnts[i]) == meshlet.vertices.end()) {
			++newvertex;
		}
	}

	// Will this triangle fit?
	if (meshlet.vertices.size() + newvertex > maxv)
		return false;

	for (uint i = 0; i < 3u; ++i) {
		if (meshlet.vertices.find(pnts[i]) == meshlet.vertices.end()) {
			meshlet.vertices.insert(pnts[i]);
		}
	}
	
	assert(meshlet.faces.find(faceid) == meshlet.faces.end());
	meshlet.faces.insert(faceid);

	return true;

}

void NewCluster::CalculateCenterRadius(const std::vector<vec3>& positions, float& radius, vec3& center)
{
	uint32_t minAxis[3] = { 0, 0, 0 };
	uint32_t maxAxis[3] = { 0, 0, 0 };

	for (uint32_t i = 1; i < positions.size(); ++i)
	{
		const vec3& pnt = positions[i];

		for (uint32_t j = 0; j < 3; ++j)
		{
			const vec3& minpnt = positions[minAxis[j]];
			const vec3& maxpnt = positions[maxAxis[j]];

			minAxis[j] = pnt[j] < minpnt[j] ? i : minAxis[j];
			maxAxis[j] = pnt[j] > maxpnt[j] ? i : maxAxis[j];
		}
	}

	float maxdis = 0.0f;
	uint axis = 0;

	for (uint i = 0; i < 3u; ++i) {
		const vec3& minpnt = positions[minAxis[i]];
		const vec3& maxpnt = positions[maxAxis[i]];

		float dis = glm::distance(minpnt, maxpnt);
		if (dis > maxdis) {
			maxdis = dis;
			axis = i;
		}
	}

	const vec3& p1 = positions[minAxis[axis]];
	const vec3& p2 = positions[maxAxis[axis]];

	vec3 tempcenter = (p1 + p2) * 0.5f;
	float tempradius = glm::length(p2 - p1) * 0.5f;
	float tempradiusSq = tempradius * tempradius;

	for (const auto& pnt : positions) {
		float dist = glm::distance(pnt, tempcenter);
		float distSq = dist * dist;

		if (distSq > tempradiusSq) {
			float k = (tempradius / dist) * 0.5f + 0.5;

			tempcenter = tempcenter * k + pnt * (1.0f - k);
			tempradius = (tempradius + dist) * 0.5f;
		}
	}

	radius = tempradius;
	center = tempcenter;
}

float NewCluster::ComputeScore(const Meshlet_built* curr, const vec3& center, const float radius, const vec3& normal, const std::array<vec3, 3>& tri, const std::array<uint, 3>& triidxs)
{
	const float reuseWeight = 0.334f;
	const float locWeight = 0.333f;
	const float oriWeight = 0.333f;
	const Meshlet_built& meshlet = *curr;
	uint reuse = 0;
	for (const auto& idx : triidxs)
		if (meshlet.vertices.find(idx) != meshlet.vertices.end())
			reuse++;
	float reuseScore = 1.0f - static_cast<float>(reuse) / 3.0f;

	float maxdis = 0.0;
	for (const auto& pnt : tri) {
		float dis = glm::distance(center, pnt);
		if (dis > maxdis)
			maxdis = dis;
	}
	float locScore = glm::log((maxdis * maxdis) / (radius * radius) + 1.0f);

	const vec3 edge1 = tri[1] - tri[0];
	const vec3 edge2 = tri[2] - tri[1];
	vec3 trinormal = glm::normalize(glm::cross(edge1, edge2));
	float oriScore = (1.0f - 1.0f * glm::dot(trinormal, normal)) / 2.0f;

	return reuseScore * reuseWeight + locWeight * locScore + oriWeight * oriScore;
}

bool NewCluster::IsMeshletFull(const Meshlet_built* curr, const int maxv, const int maxf)
{
	const Meshlet_built& meshlet = *curr;
	assert(meshlet.vertices.size() <= maxv);
	assert(meshlet.faces.size() <= maxf);

	return (meshlet.vertices.size() == maxv)
		|| (meshlet.faces.size() == maxf);
}

bool CompareScores(const std::pair<uint32_t, float>& a, const std::pair<uint32_t, float>& b)
{
	return a.second > b.second;
}

void NewCluster::PackOldmeshlets()
{
	oldmeshlets.resize(mymeshlets.size());
	for (int i = 0; i < mymeshlets.size(); ++i) {
		for (auto& face : mymeshlets[i].faces)
			oldmeshlets[i].push_back(face);
	}
}

void NewCluster::GenerateEwire(const MyMesh& mesh) {
	//我叼你大爷，出现了环的情况
	std::vector<uint> facemeshletid;
	facemeshletid.resize(nfaces);
	std::vector<uint> vertexmeshletcnt;
	vertexmeshletcnt.resize(nvertices);


	for (uint i = 0; i < oldmeshlets.size(); i++)
		for (auto& face : oldmeshlets[i])
			facemeshletid[face] = i;

	for(const auto&ver:mesh.vertices()){
		uint vidx = ver.idx();
		vertexmeshletcnt[vidx] = CalculateVertexCnt(mesh, vidx, facemeshletid);
	}

	//

	buildsets.resize(mymeshlets.size());

	for (uint i = 0; i < mymeshlets.size(); ++i) {
		//处理单个meshlet的所有情况
		for (const auto& ver : mymeshlets[i].vertices) {
			if (vertexmeshletcnt[ver] > 2)
				buildsets[i].cornerset.insert(ver);
			if (vertexmeshletcnt[ver] >= 2)
				buildsets[i].borderset.insert(ver);
		}
		//被包围的
		if (buildsets[i].cornerset.empty()) {
			buildsets[i].sign = Surrouned;
			lateaddress_surred.push_back(i);
		}

		if (buildsets[i].cornerset.size() == 1) {
			buildsets[i].sign = Surrouned;
			lateaddress_surred.push_back(i);
			std::cout << "strange topo occur" << std::endl;
		}
	}


	//由两个面判断的bord edge set总是对的
	for (uint i = 0; i < mymeshlets.size(); ++i) 
		for (const auto& face : mymeshlets[i].faces) {
			auto fh = mesh.face_handle(face);
			for (auto fhe_it = mesh.cfh_iter(fh); fhe_it.is_valid(); ++fhe_it) {
				auto opp_fidx = fhe_it->opp().face().idx();
				if (opp_fidx < 0)
					continue;
				if (mymeshlets[i].faces.find(opp_fidx) == mymeshlets[i].faces.end())
					buildsets[i].halfedgeset.insert(fhe_it->idx());
			}
		}
	
	for (auto& meshletid : lateaddress_surred) {
		auto bordedgeh = mesh.halfedge_handle(*buildsets[meshletid].halfedgeset.begin());
		auto oppfaceidx = mesh.opposite_face_handle(bordedgeh).idx();
		auto adjmeshletid = facemeshletid[oppfaceidx];
		buildsets[adjmeshletid].sign = Include;

		lateaddress_inc.push_back(adjmeshletid);

		uint vertexseed = *buildsets[meshletid].borderset.begin();
		buildsets[meshletid].cornerset.insert(vertexseed);
		buildsets[adjmeshletid].cornerset.insert(vertexseed);
	}


	//meshlet id 为25,22的时候出现了问题，可能是生成meshlet时候出现的问题
	int debugidx = 0;
	for (auto& buildset : buildsets) {
		if(buildset.sign == Normal)
			BuildEwiresOfSingleMeshlet(buildset, mesh);
		
		if (buildset.sign == Surrouned) {

			std::vector<uint> halfedgecopy;
			for (auto& edge : buildset.halfedgeset) {
				auto oppheidx = mesh.opposite_halfedge_handle(mesh.halfedge_handle(edge)).idx();
				halfedgecopy.push_back(oppheidx);
			}

			auto heh = mesh.halfedge_handle(halfedgecopy.front());
			auto oppfaceidx = mesh.face_handle(heh).idx();
			auto adjmeshletid = facemeshletid[oppfaceidx];

			assert(buildset.cornerset.size() == 1);
			BuildEwiresOfSingleMeshlet(buildset, mesh);
			assert(buildset.ewires.size() == 1);

			for (auto& he : halfedgecopy) {
				auto num_erase = buildsets[adjmeshletid].halfedgeset.erase(he);
				assert(num_erase > 0);
			}
			std::vector<uint> onewire;
			for (auto rit = buildset.ewires[0].rbegin(); rit != buildset.ewires[0].rend(); ++rit)
				onewire.push_back(*rit);
			buildsets[adjmeshletid].ewires.push_back(onewire);
		}

		if (buildset.sign == Include) {
		// do nothing
		}
		++debugidx;
	}

	for (auto& setid : lateaddress_inc) {
		auto& buildset = buildsets[setid];
		assert(buildset.sign == Include);
		assert(buildset.ewires.size() > 0);

		BuildEwiresOfSingleMeshlet(buildset, mesh);
	}

	//建立global ewire 和 meshlet ewireid map

	BuildGloEwires(buildsets);


}

uint NewCluster::CalculateVertexCnt(const MyMesh& mesh, uint vidx, const std::vector<uint>& facemeshletid)
{
	auto vh = mesh.vertex_handle(vidx);
	std::set<uint> meshletset;
	for (MyMesh::ConstVertexFaceIter vf_it = mesh.cvf_iter(vh); vf_it.is_valid(); ++vf_it)
	{
		auto fidx = vf_it->idx();
		meshletset.insert(facemeshletid[fidx]);
	}

	return meshletset.size();
}

void NewCluster::BuildEwiresOfSingleMeshlet(EwireBuildSet& buildset, const MyMesh& mesh)
{
	std::unordered_set<uint> edges = buildset.halfedgeset;
	const std::unordered_set<uint>& corners = buildset.cornerset;

	while (!edges.empty())
	{
		std::deque<uint> edgewire;
		uint startidx;
		if(buildset.ewires.empty())
			startidx = *edges.begin();
		else {
			startidx = uint(-1);
			auto lastvidx = buildset.ewires.back().back();
			auto vh = mesh.vertex_handle(lastvidx);
			for (auto voh_it = mesh.cvoh_iter(vh); voh_it.is_valid(); ++voh_it) 
				if (edges.find(voh_it->idx()) != edges.end()) {
					startidx = voh_it->idx();
					break;
				}
			if (startidx == uint(-1)) {
				//std::cout << "error may occur in build ewires" << std::endl;
				//cluster是一个环的情况
				startidx = *edges.begin();
			}
		}
		edges.erase(startidx);
		edgewire.push_back(startidx);
		//push back front可以优化是否要遍历opp
		PushEdgeAtback(edgewire, mesh, edges, corners);
		PushEdgeAtFront(edgewire, mesh, edges, corners);
		PackAndCheck(buildset.ewires, corners, edgewire, mesh);
	}

	//ReorderAndCheck(buildset);

}


void NewCluster::PushEdgeAtback(std::deque<uint>& edgewire, const MyMesh& mesh, std::unordered_set<uint>& edges, const std::unordered_set<uint> corners)
{
	do {
		auto heidx = edgewire.back();
		auto heh = mesh.halfedge_handle(heidx);
		auto toverh = mesh.to_vertex_handle(heh);
		auto vidx = toverh.idx();

		if (corners.find(vidx) != corners.end())return;

		for (MyMesh::ConstVertexOHalfedgeIter voh_it(mesh, toverh); voh_it.is_valid(); ++voh_it)
		{
			MyMesh::HalfedgeHandle temp = *voh_it;

			if (edges.find(temp.idx()) != edges.end()) {
				edges.erase(temp.idx());
				edgewire.push_back(temp.idx());
				break;
			}
			//int oppidx = voh_it->opp().idx();
			//if (edges.find(oppidx) != edges.end()) {
			//	edges.erase(oppidx);
			//	edgewire.push_back(oppidx);
			//	break;
			//}

		}

	} while (!edges.empty());

}

void NewCluster::PushEdgeAtFront(std::deque<uint>& edgewire, const MyMesh& mesh, std::unordered_set<uint>& edges, const std::unordered_set<uint> corners)
{
	do {
		auto heidx = edgewire.front();
		auto heh = mesh.halfedge_handle(heidx);
		auto fromverh = mesh.from_vertex_handle(heh);
		auto vidx = fromverh.idx();

		if (corners.find(vidx) != corners.end())return;

		for (MyMesh::ConstVertexIHalfedgeIter vih_it(mesh, fromverh); vih_it.is_valid(); ++vih_it)
		{
			MyMesh::HalfedgeHandle temp = *vih_it;

			if (edges.find(temp.idx()) != edges.end()) {
				edges.erase(temp.idx());
				edgewire.push_front(temp.idx());
				break;
			}
			//int oppidx = vih_it->opp().idx();
			//if (edges.find(oppidx) != edges.end()) {
			//	edges.erase(oppidx);
			//	edgewire.push_front(oppidx);
			//	break;
			//}
		}

	} while (!edges.empty());

}

void NewCluster::PackAndCheck(std::vector<std::vector<uint>>& ewires, const std::unordered_set<uint>& corners, std::deque<uint>& edgewire, const MyMesh& mesh)
{
	auto startedge = edgewire.front();
	auto endedge = edgewire.back();
	auto starth = mesh.halfedge_handle(startedge);
	auto endh = mesh.halfedge_handle(endedge);
	auto startvidx = mesh.from_vertex_handle(starth).idx();
	auto endvidx = mesh.to_vertex_handle(endh).idx();

	if (corners.find(startvidx) == corners.end())std::cout << "error in check pack and check";
	if (corners.find(endvidx) == corners.end())std::cout << "error in check pack and check";

	std::vector<uint> res;
	res.push_back(startvidx);

	for (auto& elem : edgewire) {
		auto edgehandle = mesh.halfedge_handle(elem);
		auto fromidx = mesh.from_vertex_handle(edgehandle).idx();
		auto toidx = mesh.to_vertex_handle(edgehandle).idx();
		if (fromidx != res.back())std::cout << "error in check pack and check";
		res.push_back(toidx);
	}

	ewires.push_back(res);
}

void NewCluster::BuildGloEwires(std::vector<EwireBuildSet>& buildsets)
{
	std::map<std::array<uint, 4>, uint> ewiremap;
	
	uint wireidx = 0;

	for (auto& buildset : buildsets) {
		buildset.ewiresid.resize(buildset.ewires.size());
		buildset.ewirereverse.resize(buildset.ewires.size());

		for(uint i=0;i<buildset.ewires.size();++i){
			auto& wire = buildset.ewires[i];
			const std::array<uint, 4> hasharray{wire[0],wire[1],wire[wire.size()-2],wire.back()};
			if (ewiremap.find(hasharray) != ewiremap.end()) {
				std::cout << "error may occur in ewire map" << std::endl;
				buildset.ewiresid[i] = ewiremap[hasharray];
				buildset.ewirereverse[i] = false;
				continue;
			}
			const std::array<uint, 4> revhasharray{hasharray[3],hasharray[2],hasharray[1],hasharray[0]};
			if (ewiremap.find(revhasharray) != ewiremap.end()) {
				buildset.ewiresid[i] = ewiremap[revhasharray];
				buildset.ewirereverse[i] = true;
				continue;
			}

			buildset.ewiresid[i] = wireidx;
			buildset.ewirereverse[i] = false;
			gloewires.push_back(wire);
			ewiremap.emplace(hasharray, wireidx++);
		}
	}

	assert(gloewires.size() == wireidx);


}

void NewCluster::ReorderAndCheck(EwireBuildSet& buildset)
{
	std::vector<std::vector<uint>> orderedwires;
	std::vector<bool> checked;
	checked.resize(buildset.ewires.size());

	orderedwires.push_back(buildset.ewires.front());
	checked[0] = true;
	while (orderedwires.size() != buildset.ewires.size()) {
		//最后一个wire的最后一个vertex的idx是下一个wire的起始idx
		uint startidx = orderedwires.back().back();
		bool flag = false;
		for(uint i=0;i<buildset.ewires.size();++i)
			if ((checked[i] == false)&&(buildset.ewires[i].front()==startidx) ){
				checked[i] = true;
				orderedwires.push_back(buildset.ewires[i]);
				flag = true;
				break;
			}
		if (flag == false) {
			std::cout << "error in reorder the ewire" << std::endl;
			break;
		}
	}
	buildset.ewires.clear();
	buildset.ewires = orderedwires;
}

bool NewCluster::SharpTriCut(const MyMesh& mesh)
{
	// mymeshlets的face和vertices都已经初始化好了

	std::vector<uint> face2meshlet;
	face2meshlet.resize(mesh.n_faces());

	for (uint mid = 0; mid < mymeshlets.size(); ++mid) {
		const auto& meshlet = mymeshlets[mid];
		for (const auto& face : meshlet.faces)
			face2meshlet[face] = mid;
	}

	struct FaceCandidate
	{
		uint faceid;
		float cost = 0;
	};

	std::vector<FaceCandidate> candidates;
	//初始化每一个triangle 对应的meshlet id
	std::vector<uint> SingleTriMeshlets;
	std::vector<uint> NeighBoor;

	std::vector<Eigen::Vector3d> centroids;
	std::vector<Eigen::Vector3d> normals;

	for (uint mid = 0; mid < mymeshlets.size(); ++mid) {
		// 计算单个meshlet的拟合平面
		Eigen::MatrixXd faceVertices(3,mymeshlets[mid].faces.size()*3);
		Eigen::Vector3d centroid;
		Eigen::Vector3d normal;

		uint vertexIdx = 0;

		for (auto& face : mymeshlets[mid].faces) {
			for (auto fvit = mesh.cfv_iter(mesh.face_handle(face)); fvit.is_valid(); ++fvit) {
				auto pnt = mesh.point(*fvit);
				faceVertices.col(vertexIdx) << pnt[0], pnt[1], pnt[2];
				++vertexIdx;
			}
		}
		fitPlane(faceVertices, centroid, normal);
		// 输出拟合平面的法向量


		centroids.push_back(centroid);
		normals.push_back(normal);

		lines.emplace_back(centroid[0], centroid[1], centroid[2]);
		lines.emplace_back(centroid[0] + normal[0]*5, centroid[1] + normal[1] * 5, centroid[2] + normal[2] * 5);
	}

	// for_each meshlet : calculate the potential triangle, store
	// for_each meshlet : for each potential triangle ,calculate the cost
	// reorder

	for (uint mid = 0; mid < mymeshlets.size(); ++mid) {
		for (auto& face : mymeshlets[mid].faces) {
			std::vector<uint> meshletsids;
			for (auto ffit = mesh.cff_iter(mesh.face_handle(face)); ffit.is_valid(); ++ffit)
				if(face2meshlet[ffit->idx()]!=mid)
					meshletsids.push_back(face2meshlet[ffit->idx()]);
			
			if ((meshletsids.size() == 2) && (meshletsids[0] == meshletsids[1])) {
				FaceCandidate temp;
				temp.faceid = face;
				temp.cost = EvaluateCost(centroids[mid],normals[mid],face2meshlet,face,mesh);
				candidates.push_back(temp);
			}	

			if ((meshletsids.size() == 3) && (meshletsids[0] == meshletsids[1]) && (meshletsids[0] == meshletsids[2])) {
				SingleTriMeshlets.push_back(mid);
				NeighBoor.push_back(meshletsids[0]);
			}
		}
	}

	if (candidates.empty())
		return false;

	std::sort(candidates.begin(), candidates.end(), [](const FaceCandidate& a, const FaceCandidate& b) -> bool
	{
		return a.cost > b.cost;
	});

	
	for (uint i = 0; i < SingleTriMeshlets.size(); ++i) {
		assert(mymeshlets[SingleTriMeshlets[i]].faces.size() == 1);
		uint faceid = *mymeshlets[SingleTriMeshlets[i]].faces.begin();
		uint neighmeshletid = NeighBoor[i];
		
		if (mymeshlets[neighmeshletid].faces.size() >= maxf)
			continue;

		mymeshlets[neighmeshletid].faces.insert(faceid);
		mymeshlets[SingleTriMeshlets[i]].faces.clear();
		mymeshlets[SingleTriMeshlets[i]].vertices.clear();
	}


	for (auto& elem : candidates) {
		auto face = elem.faceid;
		uint mid = face2meshlet[face];
		std::vector<uint> meshletsids;
		for (auto ffit = mesh.cff_iter(mesh.face_handle(face)); ffit.is_valid(); ++ffit)
			if (face2meshlet[ffit->idx()] != mid)
				meshletsids.push_back(face2meshlet[ffit->idx()]);

		if ((meshletsids.size() == 2) && (meshletsids[0] == meshletsids[1])&& (mymeshlets[meshletsids[0]].faces.size()<maxf)) {
			face2meshlet[face] = meshletsids[0];
			mymeshlets[mid].faces.erase(face);
			mymeshlets[meshletsids[0]].faces.insert(face);
		}
		else
			continue;

		if (mymeshlets[mid].faces.empty())
			mymeshlets[mid].vertices.clear();

		for (auto fhit = mesh.cfh_iter(mesh.face_handle(face)); fhit.is_valid(); ++fhit) {
			if (face2meshlet[fhit->opp().face().idx()] == mid)
				mymeshlets[mid].vertices.erase(fhit->next().to().idx());
		}
	}

	for (auto it = mymeshlets.begin(); it != mymeshlets.end();) {
		if (it->faces.empty())
			it = mymeshlets.erase(it);
		else
			it++;
	}

	return true;
	
}

void NewCluster::fitPlane(const Eigen::MatrixXd& points, Eigen::Vector3d& centroid, Eigen::Vector3d& normal)
{
	int numpoints = points.cols();

	// 计算中心点
	centroid = points.rowwise().mean();

	// 计算协方差矩阵
	Eigen::MatrixXd centered = points.colwise() - centroid;
	Eigen::Matrix3d covariance = (centered * centered.transpose()) / static_cast<double>(numpoints);

	// 计算特征值和特征向量
	Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver(covariance);
	Eigen::Vector3d eigenvalues = solver.eigenvalues();
	Eigen::Matrix3d eigenvectors = solver.eigenvectors();

	// 选择最小特征值对应的特征向量作为法线向量
	int minIndex;
	eigenvalues.minCoeff(&minIndex);
	normal = eigenvectors.col(minIndex).normalized();
}

void NewCluster::TryShapeHeal(const MyMesh& mesh)
{
	std::vector<uint> errormeshlets;

	detecters.resize(mymeshlets.size());
	for (uint mid = 0; mid < mymeshlets.size(); ++mid) {
		bool res = ShapeErrorDetect(mesh, mid);
		if (res == false)
			errormeshlets.push_back(mid);
	}

	if (errormeshlets.size() == 0)
		return;

	do {
		for (auto& mid : errormeshlets) {
			if (detecters[mid].topotype == Torus)
				TorusHeal(mesh, mid);
			else if (detecters[mid].topotype == SharpConnect)
				NewSharpConnectHeal(mesh, mid);
			else
				std::cout << "what happened?" << std::endl;
		}

		errormeshlets.clear();

		for (uint mid = 0; mid < mymeshlets.size(); ++mid) {
			bool res = ShapeErrorDetect(mesh, mid);
			if (res == false)
				errormeshlets.push_back(mid);
		}

	} while (errormeshlets.size() != 0);

}

//如果好的返回true
bool NewCluster::ShapeErrorDetect(const MyMesh& mesh, uint mid)
{
	if (mid == 1638)
		std::cout << "debug" << std::endl;


	auto& detecter = detecters[mid];
	if (detecter.topotype != Unchecked)
		return detecter.topotype == Correct;

	detecter.halfedgeset.clear();
	detecter.specialpnts.clear();
	detecter.newfacestart.clear();

	for (auto& face : mymeshlets[mid].faces) {
		auto fh = mesh.face_handle(face);
		for (auto fhit = mesh.cfh_iter(fh); fhit.is_valid(); ++fhit) {
			auto oppfid = fhit->opp().face().idx();
			if (mymeshlets[mid].faces.find(oppfid) == mymeshlets[mid].faces.end())
				detecter.halfedgeset.insert(fhit->idx());
		}
	}

	auto& edgeset = detecter.halfedgeset;
	std::map<uint, int> vertexcnt;
	int loopcnt = 0;
	std::vector<std::vector<uint>> loops;

	//find first loop
	while (!edgeset.empty())
	{
		std::vector<uint> oneloop;
		uint seedhe = *edgeset.begin();
		edgeset.erase(seedhe);
		auto heh = mesh.halfedge_handle(seedhe);
		auto startidx = mesh.from_vertex_handle(heh).idx();
		auto nextidx = mesh.to_vertex_handle(heh).idx();
		vertexcnt[startidx] += 1;
		vertexcnt[nextidx] += 1;

		oneloop.push_back(heh.idx());

		do {
			auto vh = mesh.vertex_handle(nextidx);
			for (auto cvoit = mesh.cvoh_iter(vh); cvoit.is_valid(); ++cvoit) {
				if (edgeset.find(cvoit->idx()) != edgeset.end()) {
					edgeset.erase(cvoit->idx());
					vertexcnt[nextidx] += 1;
					if(nextidx == 1421)
						std::cout << "debug" << std::endl;
					nextidx = cvoit->to().idx();
					oneloop.push_back(cvoit->idx());
					vertexcnt[nextidx] += 1;
					break;
				}
			}
		} while (nextidx != startidx);
		loopcnt++;
		loops.push_back(oneloop);
	}

	if (loopcnt == 1)
		detecter.topotype = Correct;
	else
		detecter.topotype = Torus;

	for (auto& elem : vertexcnt)
		if (elem.second != 2)
			detecter.specialpnts.push_back(elem.first);
	if (!detecter.specialpnts.empty())
		detecter.topotype = SharpConnect;

	if (detecter.topotype == SharpConnect)
		std::cout << "Check type error mid: sharp connect   " << mid << std::endl;
	if(detecter.topotype == Torus)
		std::cout << "Check type error mid: torus           " << mid << std::endl;

	if (detecter.topotype == Torus) {
		GenerateTorusNewFaceStart(mesh, loops, mid);
	}

	return detecter.topotype==Correct;
}

void NewCluster::TorusHeal(const MyMesh& mesh, uint mid)
{
	//对mymeshlets[mid]进行重构，重新划分face,更新vertex
	auto& meshlet = mymeshlets[mid];

	auto& faceset = mymeshlets[mid].faces;
	uint totalfacenum = faceset.size();
	std::unordered_set<uint> newfaceset;

	//GrowFaceSet(mesh, newfaceset, faceset, detecters[mid].newfacestart, totalfacenum / 2);
	for (auto& face : detecters[mid].newfacestart)
		newfaceset.insert(face);
	for (auto& face : newfaceset)
		meshlet.faces.erase(face);

	//更新老的
	meshlet.vertices.clear();
	for (auto& face : meshlet.faces) {
		auto fh = mesh.face_handle(face);
		for (auto fvit = mesh.cfv_iter(fh); fvit.is_valid(); ++fvit)
			meshlet.vertices.insert(fvit->idx());
	}
	detecters[mid].topotype = Unchecked;
	detecters[mid].halfedgeset.clear();
	detecters[mid].specialpnts.clear();
	detecters[mid].newfacestart.clear();

	//更新新的
	Meshlet_built newmeshlet;
	newmeshlet.faces = newfaceset;

	for (auto& face : newmeshlet.faces) {
		auto fh = mesh.face_handle(face);
		for (auto fvit = mesh.cfv_iter(fh); fvit.is_valid(); ++fvit)
			newmeshlet.vertices.insert(fvit->idx());
	}
	ShapeHeal newdetect;
	newdetect.topotype = Unchecked;

	mymeshlets.push_back(newmeshlet);
	detecters.push_back(newdetect);

}

void NewCluster::SharpConnectHeal(const MyMesh& mesh, uint mid)
{
	//对mymeshlets[mid]进行重构，重新划分face,更新vertex
	auto& meshlet = mymeshlets[mid];
	uint startv = detecters[mid].specialpnts[0];
	auto vh = mesh.vertex_handle(startv);
	uint startf;
	auto& faceset = meshlet.faces;

	for(auto cvfit = mesh.cvf_iter(vh);cvfit.is_valid();++cvfit)
		if (faceset.find(cvfit->idx()) != faceset.end()) {
			startf = cvfit->idx();
			break;
		}

	uint totalfacenum = faceset.size();
	std::unordered_set<uint> newfaceset;

	GrowFaceSet(mesh, newfaceset, faceset, startf,totalfacenum/2);

	//更新老的
	meshlet.vertices.clear();
	for (auto& face : meshlet.faces) {
		auto fh = mesh.face_handle(face);
		for (auto fvit = mesh.cfv_iter(fh); fvit.is_valid(); ++fvit)
			meshlet.vertices.insert(fvit->idx());
	}
	detecters[mid].topotype = Unchecked;
	detecters[mid].halfedgeset.clear();
	detecters[mid].specialpnts.clear();
	detecters[mid].newfacestart.clear();

	//更新新的
	Meshlet_built newmeshlet;
	newmeshlet.faces = newfaceset;

	for (auto& face : newmeshlet.faces) {
		auto fh = mesh.face_handle(face);
		for (auto fvit = mesh.cfv_iter(fh); fvit.is_valid(); ++fvit)
			newmeshlet.vertices.insert(fvit->idx());
	}
	ShapeHeal newdetect;
	newdetect.topotype = Unchecked;

	mymeshlets.push_back(newmeshlet);
	detecters.push_back(newdetect);

}

void NewCluster::NewSharpConnectHeal(const MyMesh& mesh, uint mid)
{
	//对mymeshlets[mid]进行重构，重新划分face,更新vertex
	auto& meshlet = mymeshlets[mid];
	uint startv = detecters[mid].specialpnts[0];
	auto vh = mesh.vertex_handle(startv);

	std::set<uint> face_corner;
	auto& faceset = meshlet.faces;

	for (auto cvfit = mesh.cvf_iter(vh); cvfit.is_valid(); ++cvfit)
		if (faceset.find(cvfit->idx()) != faceset.end()) {
			face_corner.insert(cvfit->idx());
		}

	std::unordered_set<uint> newfaceset;
	uint startf = *face_corner.begin();
	newfaceset.insert(startf);
	face_corner.erase(startf);

	std::deque<uint> que;
	que.push_back(startf);

	while (!que.empty())
	{
		uint faceid = que.front();
		que.pop_front();
		auto fh = mesh.face_handle(faceid);
		for (auto ffit = mesh.cff_iter(fh); ffit.is_valid(); ++ffit) {
			if (face_corner.find(ffit->idx()) != face_corner.end()) {
				newfaceset.insert(ffit->idx());
				face_corner.erase(ffit->idx());
				que.push_back(ffit->idx());
			}
		}
	}
	
	for (auto face : newfaceset)
		meshlet.faces.erase(face);

	//更新老的
	meshlet.vertices.clear();
	for (auto& face : meshlet.faces) {
		auto fh = mesh.face_handle(face);
		for (auto fvit = mesh.cfv_iter(fh); fvit.is_valid(); ++fvit)
			meshlet.vertices.insert(fvit->idx());
	}
	detecters[mid].topotype = Unchecked;
	detecters[mid].halfedgeset.clear();
	detecters[mid].specialpnts.clear();
	detecters[mid].newfacestart.clear();

	//更新新的
	Meshlet_built newmeshlet;
	newmeshlet.faces = newfaceset;

	for (auto& face : newmeshlet.faces) {
		auto fh = mesh.face_handle(face);
		for (auto fvit = mesh.cfv_iter(fh); fvit.is_valid(); ++fvit)
			newmeshlet.vertices.insert(fvit->idx());
	}
	ShapeHeal newdetect;
	newdetect.topotype = Unchecked;

	mymeshlets.push_back(newmeshlet);
	detecters.push_back(newdetect);
}

void NewCluster::GrowFaceSet(const MyMesh& mesh, std::unordered_set<uint>& newfaceset, std::unordered_set<uint>& faceset, uint seedface, uint cutlimit)
{
	std::deque<uint> fifo;
	fifo.push_back(seedface);

	while (newfaceset.size() < cutlimit) {
		uint faceid = fifo.front();
		fifo.pop_front();

		if (newfaceset.find(faceid) != newfaceset.end())
			continue;

		faceset.erase(faceid);
		newfaceset.insert(faceid);


		auto fh = mesh.face_handle(faceid);
		for (auto ffit = mesh.cff_iter(fh); ffit.is_valid(); ++ffit) {
			if (faceset.find(ffit->idx()) != faceset.end())
				fifo.push_back(ffit->idx());
		}
	}

}

void NewCluster::GrowFaceSet(const MyMesh& mesh, std::unordered_set<uint>& newfaceset, std::unordered_set<uint>& faceset, std::vector<uint> seedfaces, uint cutlimit)
{
	std::deque<uint> fifo;
	for(auto&face:seedfaces)
		fifo.push_back(face);

	while (newfaceset.size() < cutlimit) {
		uint faceid = fifo.front();
		fifo.pop_front();

		if (newfaceset.find(faceid) != newfaceset.end())
			continue;

		faceset.erase(faceid);
		newfaceset.insert(faceid);


		auto fh = mesh.face_handle(faceid);
		for (auto ffit = mesh.cff_iter(fh); ffit.is_valid(); ++ffit) {
			if (faceset.find(ffit->idx()) != faceset.end())
				fifo.push_back(ffit->idx());
		}
	}
}

void NewCluster::GenerateTorusNewFaceStart(const MyMesh& mesh, const std::vector<std::vector<uint>>& loops, uint mid)
{
	auto& detector = detecters[mid];
	auto& meshlet = mymeshlets[mid];

	uint he = loops[0][0];
	auto heh = mesh.halfedge_handle(he);
	auto startfh = mesh.face_handle(heh);
	uint startidx = startfh.idx();

	std::unordered_set<uint> targetfaces;
	for (auto& heid : loops[1]) {
		auto heidh = mesh.halfedge_handle(heid);
		auto fh = mesh.face_handle(heidh);
		targetfaces.insert(fh.idx());
	}

	std::vector<std::vector<std::pair<uint,int>>> depthrecord;


	auto res = BFS(mesh, startidx, targetfaces,mymeshlets[mid].faces);
	for (auto& face : res)
		detector.newfacestart.push_back(face);
}

std::vector<uint> NewCluster::BFS(const MyMesh& mesh, uint startidx, const std::unordered_set<uint>& targets, const std::unordered_set<uint>& faceset)
{
	std::unordered_set<uint> visited;

	std::vector<std::vector<std::pair<uint, int>>> depthrecord;
	std::vector<std::pair<uint, int>> startdepth{ std::pair<uint, int>{startidx, -1}};
	depthrecord.push_back(startdepth);
	visited.insert(startidx);
	bool sign = false;
	uint mytarget;
	uint layer = 0;
	uint pre;

	do {
		auto& currentdepth = depthrecord[layer];
		std::vector<std::pair<uint, int>> nextdepth;

		for (auto& elem : currentdepth) {
			auto fidx = elem.first;
			auto fh = mesh.face_handle(fidx);
			for (auto ffit = mesh.cff_iter(fh); ffit.is_valid(); ++ffit) {
				if (faceset.find(ffit->idx()) == faceset.end())
					continue;

				if (visited.find(ffit->idx()) == visited.end()) {
					nextdepth.push_back(std::pair<uint, int>{ffit->idx(), elem.first});
					visited.insert(ffit->idx());
				}
				if (targets.find(ffit->idx()) != targets.end()) {
					sign = true;
					mytarget = ffit->idx();

					break;
				}
			}
			if (sign == true)
				break;
		}
		depthrecord.push_back(nextdepth);
		layer++;

	} while (sign == false);
	std::vector<uint> res{ mytarget };

	for (uint i = depthrecord.size() - 1; i != 0; --i) {
		for(auto&elem:depthrecord[i])
			if (elem.first == res.back()) {
				res.push_back(elem.second);
				break;
			}
	}
	assert(res.back() == startidx);
	return res;

}

double NewCluster::point2PlaneDist(const Eigen::Vector3d& point, const Eigen::Vector3d& centroid, const Eigen::Vector3d& normal)
{
	Eigen::Vector3d diff = point - centroid;
	double dist = std::abs(diff.dot(normal));
	return dist;
}

double NewCluster::EvaluateCost(const Eigen::Vector3d& centroid, const Eigen::Vector3d& normal, const std::vector<uint>& meshletid, uint faceid, const MyMesh& mesh)
{
	double cost = 0;
	auto fh = mesh.face_handle(faceid);
	for (auto cfvit = mesh.cfv_iter(fh); cfvit.is_valid(); ++cfvit) {
		auto pnt = mesh.point(*cfvit);
		cost += point2PlaneDist(Eigen::Vector3d{ pnt[0],pnt[1],pnt[2] }, centroid, normal);
	}

	for (auto cfheit = mesh.cfh_iter(fh); cfheit.is_valid(); ++cfheit) {
		auto oppfaceid = cfheit->opp().face().idx();
		auto pnt0 = mesh.point(cfheit->from());
		auto pnt1 = mesh.point(cfheit->to());
		vec3 p0{ pnt0[0],pnt0[1],pnt0[2] };
		vec3 p1{ pnt1[0],pnt1[1],pnt1[2] };
		double dis = glm::distance(p0, p1);

		if (meshletid[oppfaceid] != meshletid[faceid])
			cost += dis;
		else
			cost -= dis;
	}
	//cost越高 该faceid的毛刺程度，偏离程度越大
	return cost;
}
