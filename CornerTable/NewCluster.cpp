#include "NewCluster.h"

bool CompareScores(const std::pair<uint32_t, float>& a, const std::pair<uint32_t, float>& b);

NewCluster::NewCluster(uint maxverts, uint maxtris, const MyMesh& mesh)
{
	nfaces = mesh.n_faces();
	nvertices = mesh.n_vertices();

	mymeshlets.clear();
	mymeshlets.emplace_back();
	Meshlet_built* curr = &mymeshlets.back();

	std::vector<std::pair<uint, float>> candidates;
	std::unordered_set<uint> candidateCheck;
	std::vector<bool> checklist(nfaces,false);
	

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
	std::vector<EwireBuildSet> buildsets;
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