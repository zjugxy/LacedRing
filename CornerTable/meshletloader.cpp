#include"Meshlet.h"
#include<random>

//TS
void MeshletLoad(std::vector<TS_meshlet>& loader, MyMesh mesh, Meshlets meshlets,std::vector<vec3>& geometryinfo) {

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
				geometryinfo.emplace_back(pnt[0], pnt[1], pnt[2]);
				++cnt;
			}
		}
		//单个meshlet遍历完成
		temp.primitive_cnt = pricnt;
		temp.color = vec3(dis(gen), dis(gen), dis(gen));
		loader.push_back(temp);
	}

}