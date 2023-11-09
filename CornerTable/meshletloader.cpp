#include"Meshlet.h"
#include<random>
#include<unordered_map>
//TS
void TS_MeshletLoad(std::vector<TS_meshlet>& loader, MyMesh mesh, Meshlets meshlets,std::vector<vec4>& geometryinfo) {

	// 使用随机设备作为种子
	std::random_device rd;
	// 使用随机设备生成引擎
	std::mt19937 gen(rd());
	// 创建均匀分布对象，范围在 0 到 1 之间
	std::uniform_real_distribution<float> dis(0.0, 1.0);

	geometryinfo.clear();
	loader.clear();
	int cnt = 0;

	for (const auto& meshlet : meshlets) {
		TS_meshlet temp;
		temp.vertex_begin = cnt;
		int pricnt = 0;
		for (const auto& faceidx : meshlet) {
			pricnt++;
			MyMesh::FaceHandle fh = mesh.face_handle(faceidx);
			for (MyMesh::FaceVertexIter fv_it = mesh.fv_iter(fh); fv_it.is_valid(); ++fv_it) {
				auto vh = *fv_it;
				auto pnt = mesh.point(vh);
				geometryinfo.emplace_back(pnt[0], pnt[1], pnt[2],1.0f);
				++cnt;
			}
		}
		//单个meshlet遍历完成
		temp.primitive_cnt = pricnt;
		temp.color = vec4(dis(gen), dis(gen), dis(gen),1.0f);
		loader.push_back(temp);
	}

}

void IX_meshletLoad(std::vector<IX_meshlet>& loader, MyMesh mesh, Meshlets meshlets, std::vector<vec4>& geoinfo, std::vector<int>& vertexidx, std::vector<int>& primidx)
{
	// 使用随机设备作为种子
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<float> dis(0.0, 1.0);

	geoinfo.clear();
	loader.clear();
	vertexidx.clear();
	primidx.clear();
	int cnt = 0;

	//load geoinfo
	for (auto& ver : mesh.vertices()) {
		auto pnt = mesh.point(ver);
		geoinfo.emplace_back(pnt[0],pnt[1],pnt[2],1.0f);
	}

	int vertex_glo_begin = 0;
	int prim_glo_begin = 0;


	for (const auto& meshlet : meshlets) {
		IX_meshlet temp;
		temp.vertex_begin = vertex_glo_begin;
		temp.primbegin = prim_glo_begin;
		//create a map vhidx --> primitiveidx
		std::unordered_map<int, int> vh2prim;
		int startidx = 0;


		for (const auto& faceidx : meshlet) {
			MyMesh::FaceHandle fh = mesh.face_handle(faceidx);
			for (MyMesh::FaceVertexIter fv_it = mesh.fv_iter(fh); fv_it.is_valid(); ++fv_it) {
				auto vh = *fv_it;
				auto vhidx = vh.idx();
				if (vh2prim.find(vhidx) == vh2prim.end()) {
					vh2prim[vhidx] = startidx++;
					vertexidx.push_back(vhidx);
				}
				primidx.push_back(vh2prim[vhidx]);
			}
		}
		temp.vertex_cnt = vertexidx.size() - vertex_glo_begin;
		vertex_glo_begin = vertexidx.size();
		temp.primcnt = (primidx.size() - prim_glo_begin)/3;
		prim_glo_begin = primidx.size();
		temp.color = vec4(dis(gen), dis(gen), dis(gen), 1.0f);
		loader.push_back(temp);

	}

}

void SC_meshletLoad(std::vector<SC_meshlet>& loader, MyMesh mesh, Meshlets meshlets, std::vector<vec4>& geoinfo, std::vector<int>& primidx)
{

	// 使用随机设备作为种子
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<float> dis(0.0, 1.0);

	geoinfo.clear();
	loader.clear();
	primidx.clear();
	int cnt = 0;

	//load geoinfo


	int vertex_glo_begin = 0;
	int prim_glo_begin = 0;


	for (const auto& meshlet : meshlets) {
		SC_meshlet temp;
		temp.vertex_begin = vertex_glo_begin;
		temp.primbegin = prim_glo_begin;
		//create a map vhidx --> primitiveidx
		std::unordered_map<int, int> vh2prim;
		int startidx = 0;
		for (const auto& faceidx : meshlet) {
			MyMesh::FaceHandle fh = mesh.face_handle(faceidx);
			for (MyMesh::FaceVertexIter fv_it = mesh.fv_iter(fh); fv_it.is_valid(); ++fv_it) {
				auto vh = *fv_it;
				auto vhidx = vh.idx();
				if (vh2prim.find(vhidx) == vh2prim.end()) {
					vh2prim[vhidx] = startidx++;
					auto pnt = mesh.point(vh);
					geoinfo.emplace_back(pnt[0], pnt[1], pnt[2], 1.0);
				}
				primidx.push_back(vh2prim[vhidx]);
			}
		}
		temp.vertex_cnt = geoinfo.size() - vertex_glo_begin;
		vertex_glo_begin = geoinfo.size();
		temp.primcnt = (primidx.size() - prim_glo_begin) / 3;
		prim_glo_begin = primidx.size();
		temp.color = vec4(dis(gen), dis(gen), dis(gen), 1.0f);
		loader.push_back(temp);
	}
}

