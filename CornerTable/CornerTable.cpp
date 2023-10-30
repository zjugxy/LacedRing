#include "CornerTable.h"
#include <OpenMesh/Core/Utils/PropertyManager.hh>
#include<unordered_map>
#include<list>

// ----------------------------------------------------------------------------

typedef OpenMesh::TriMesh_ArrayKernelT<>  MyMesh;


CornerTable::CornerTable(MyMesh input)
{
	//corner数量等于3倍的face
	Vtable = std::vector<int>(input.n_faces() * 3, -1);
	Otable = std::vector<int>(input.n_faces() * 3, -1);
	Ctable = std::vector<int>(input.n_vertices(), -1);

	for (const auto& ver : input.vertices()) {
		auto pnt = input.point(ver);
		vertices.push_back(vec3{pnt[0], pnt[1], pnt[2]});
		if (ver.is_boundary()) {
			std::cout << "the model is not a 2-manifold" << std::endl;
			return;
		}
	}

	//初始化每条halfedge对应的corner，并且填充Vtable
	auto cornerindex = OpenMesh::HProp<int>(input);
	for (const auto& face : input.faces()) {
		int start = face.idx() * 3+2;
		for (const auto& he : face.halfedges()) {
			cornerindex[he] = start;
			Vtable[start] = he.next().to().idx();
			start--;
		}
	}

	//初始化Otable
	for (const auto& face : input.faces()) {
		for (const auto& he : face.halfedges()) {
			int corneridx = cornerindex[he];
			auto opphe = he.opp();
			Otable[corneridx] = cornerindex[opphe];
		}
	}

	//初始化Ctable
	for (int i = 0; i < Vtable.size(); i++) {
		Ctable[Vtable[i]] = i;
	}

	//Otable check
	for (const auto& face : input.faces()) {
		//std::cout << "face id" << face.idx()<<std::endl;
		for (const auto& he : face.halfedges()) {
			int corneridx = cornerindex[he];
			//std::cout << "corner index is " << corneridx << " its opp is " << Otable[corneridx] << " its vertex is " << Vtable[corneridx] << std::endl;
			if (corneridx != Otable[Otable[corneridx]]) {
				std::cout << "error" << std::endl;
				break;
			}
		}
	}

	//Vtable check using checksum
	for (const auto& face : input.faces()) {
		int checksum = 0;
		for (const auto& ver : face.vertices())checksum += ver.idx();
		for (const auto& he : face.halfedges()) {
			int temp = cornerindex[he];
			checksum -= Vtable[temp];
		}
		if (checksum)
		{
			std::cout << "error" << std::endl;
			break;
		}
	}

	////ct information print
	//for (const auto& face : input.faces()) {
	//	std::cout << "face id" << face.idx()<<std::endl;
	//	for (const auto& he : face.halfedges()) {
	//		int corneridx = cornerindex[he];
	//		std::cout << "corner index: " << corneridx << " opp: " << Otable[corneridx] << " vertex: " << Vtable[corneridx] << std::endl;
	//	}
	//}

}

vec3 CornerTable::getVertexValue(int corner)const
{
	return vertices[Vertex(corner)];
}

int CornerTable::Vertex(int corner)const
{
	return Vtable[corner];
}

int CornerTable::Tri(int corner)const
{
	return corner/3;
}

int CornerTable::Next(int corner)const
{
	int base = 3 * (corner / 3);
	return base+(corner+1)%3;
}

int CornerTable::Prev(int corner)const
{
	int base = 3 * (corner / 3);
	return base + (corner + 2) % 3;
}

int CornerTable::Left(int corner)const
{
	return Opp(Next(corner));
}

int CornerTable::Right(int corner)const
{
	return Opp(Prev(corner));
}

int CornerTable::Opp(int corner)const
{
	return Otable[corner];
}

int CornerTable::CornerofVetex(int vertexidx)const
{
	return Ctable[vertexidx];
}

int CornerTable::Swin(int corner) const
{
	return Next(Left(corner));
}

void CornerTable::HamiltonianCycle(int start)
{
	Mvertex = std::vector<int>(vertices.size(), 0);
	Mtri = std::vector<int>(Vtable.size() / 3, 0);

	int corner = start;
	int end = Opp(start);
	Mvertex[Vertex(Prev(corner))] = 1;
	Mvertex[Vertex(Next(corner))] = 1;
	int vercnt = 2;
	int tricnt = 0;

	std::vector<int> listcorner1;
	std::vector<int> listcorner2;
	listcorner1.push_back(start);
	listcorner2.push_back(Prev(start));
	listcorner2.push_back(Prev(Prev(start)));


	do {
		if (Mvertex[Vertex(corner)] == 0) {
			Mvertex[Vertex(corner)] = Mtri[Tri(corner)] = 1;
			vercnt++;
			tricnt++;
			listcorner1.push_back(corner);

		}
		else if (Mtri[Tri(corner)] == 0) {
			listcorner2.push_back(listcorner1.back());
			listcorner1.pop_back();
			corner = Opp(corner);
			
			//std::cout <<"O  " << Vertex(corner) << std::endl;

		}
		corner = Right(corner);
		
		//std::cout <<"R  " << Vertex(corner) << std::endl;
	
	} while (corner != end);

	listcorner2.pop_back();
	listcorner2.push_back(listcorner2.front());

	//for (auto& elem : listcorner1)
	//	std::cout << elem << std::endl;
	//std::cout << "##############" << std::endl;

	//for (auto& elem : listcorner2)
	//	std::cout << elem << std::endl;

	ringcorner = listcorner2;
	for (auto& elem : ringcorner)
		ringvertex.push_back(Vertex(elem));

	std::cout << "ring vertex num" << ringvertex.size() - 1 << std::endl;//最后一个加上了初始点
	std::cout << "tri marked: " << tricnt<<std::endl;
	std::cout << "ver marked: " << vercnt << std::endl;

	if (vertices.size() != ringvertex.size() - 1)
		std::cout << "Error in not every point in Ring!!!" << std::endl;

	TriangleClassification();
}




void CornerTable::BetterHamiltonianCycle(int start) {
	Mvertex = std::vector<int>(vertices.size(), 0);
	Mtri = std::vector<int>(Vtable.size() / 3, 0);

	int corner = start;
	int end = Opp(start);
	Mvertex[Vertex(Prev(corner))] = 1;
	Mvertex[Vertex(Next(corner))] = 1;
	int vercnt = 2;
	int tricnt = 0;

	std::list<int> vertexidxlist;
	vertexidxlist.push_back(Vertex(Next(start)));
	vertexidxlist.push_back(Vertex(Prev(start)));
	vertexidxlist.push_back(Vertex(Next(start)));
	std::list<int>::iterator ptr;
	ptr = vertexidxlist.begin();
	std::advance(ptr, 2);

	do {
		if (Mvertex[Vertex(corner)] == 0) {
			Mvertex[Vertex(corner)] = Mtri[Tri(corner)] = 1;
			vercnt++;
			tricnt++;
			//三角形对应三个顶点，其中有一个边（两个顶点已经在vertexidxlist中）
			//找到该边，并在其中插入点
			// B,A,it,B -> B,A,C,B - > B,A,C,it,B -> B,A,C,D,B

			ptr = vertexidxlist.insert(ptr, Vertex(corner));
			//check ptr两边是不是对应的添加的三角形的两个顶点

			if (simplecheck(Tri(corner), ptr) == false) {
				std::cout << "error in Hamilton" << std::endl;
					return;
			}
			ptr++;
		}
		else if (Mtri[Tri(corner)] == 0) {
			ptr--;
			//B,A,C,it,B -> B,A,it,C,B
			//CB无法invade，尝试invade AC
			corner = Opp(corner);
		}
		corner = Right(corner);
	} while (corner != end);

	for (auto& elem : vertexidxlist)
		ringvertex.push_back(elem);
	ringvertex.pop_back();

	std::cout << "ring num = " << ringvertex.size() << std::endl;
	std::cout << "vertex marked = " << vercnt << std::endl;
}

bool CornerTable::simplecheck(int triid,  std::list<int>::iterator& it)const {
	int sum = Vertex(triid * 3) + Vertex(triid * 3 + 1) + Vertex(triid * 3 + 2);
	auto temp = it;
	int cmpsum = *it;
	temp--;
	cmpsum += *temp;
	temp++; temp++;
	cmpsum += *temp;
	return sum == cmpsum;
}


void CornerTable::TriangleClassification() {
	facerightleftwithcolor = std::vector<int>(Vtable.size() / 3, 0);

	triangletypes = std::vector<TypeTri>(Vtable.size() / 3, TypeTri::T0);
	//遍历
	int totalnum = ringvertex.size();
	for (int i = 0; i < ringvertex.size(); i++) {
		int thirdcorner = -1;
		int triidx = Vertex2Tris(ringvertex[i], ringvertex[(i + 1)%totalnum],thirdcorner);

		//v0 v1 thirdcorner 该三角形位于边的右边
		if (Next(thirdcorner) == ringvertex[i]) {
			RingRightCorners.push_back(thirdcorner);
			RingLeftCorners.push_back(Opp(thirdcorner));
			facerightleftwithcolor[thirdcorner / 3] = 1;//right 1
			facerightleftwithcolor[Opp(thirdcorner) / 3] = 2;//right 1
		}
		else {
			RingLeftCorners.push_back(thirdcorner);
			RingRightCorners.push_back(Opp(thirdcorner));
			facerightleftwithcolor[thirdcorner / 3] = 2;
			facerightleftwithcolor[Opp(thirdcorner) / 3] = 1;//right 1

		}
		Tri_Ring_Add(Tri(thirdcorner));
		Tri_Ring_Add(Tri(Opp(thirdcorner)));
	}
	//处理
	//需要遍历所有的T0三角形标注T1_ir T2_ir
	for (int i = 0; i < triangletypes.size(); i++) {
		if (triangletypes[i] == TypeTri::T0) {
			int corner0 = i * 3, corner1 = i * 3 + 1, corner2 = i * 3 + 2;
			int tri0 = Tri(Opp(corner0)), tri1 = Tri(Opp(corner1)), tri2 = Tri(Opp(corner2));
			AdjustTritype_AdjecentT0(tri0);
			AdjustTritype_AdjecentT0(tri1);
			AdjustTritype_AdjecentT0(tri2);
			T0trilist.push_back(i);
		}
	}
	//std::cout << sizeof(TypeTri) << std::endl;

	for (int i = 0; i < triangletypes.size(); i++) {
		switch (triangletypes[i])
		{
		case TypeTri::T0:
			RenderTriType.push_back(0); break;
		case TypeTri::T1:
			RenderTriType.push_back(1); break;
		case TypeTri::T2:
			RenderTriType.push_back(2); break;
		case TypeTri::T1_irregular:
			RenderTriType.push_back(3); break;
		case TypeTri::T2_irregular:
			RenderTriType.push_back(4); break;
		default:
			break;
		}
	}

}

int CornerTable::Vertex2Tris(int v0, int v1,int& corner)const {
	int c0 = Ctable[v0];
	int c1 = Ctable[v1];
	int temp = c0;

	do {
		int tri0 = Tri(temp);
		//若该三角形被访问过
		if (Mtri[tri0]) {
			//遍历v1对应的所有corner找一致的三角形
			int temp1 = c1;
			do {
				if (Tri(temp1) == tri0) {
					corner = 3 * tri0 + (3 - temp % 3 - temp1 % 3);
					return tri0;
				}
				temp1 = Swin(temp1);
			} while (temp1 != c1);
			temp = Swin(temp);
		}
		else {
			temp = Swin(temp);
		}
	} while (temp != c0);
	std::cout << "error in Vertex2Tris " << std::endl;
}

void CornerTable::Tri_Ring_Add(int triid) {
	if (triangletypes[triid] == TypeTri::T0)
		triangletypes[triid] = TypeTri::T1;
	else if (triangletypes[triid] == TypeTri::T1)
		triangletypes[triid] = TypeTri::T2;
	else {
		std::cout << "error in Tri ring add in triangle classification" << std::endl;
	}
}

void CornerTable::AdjustTritype_AdjecentT0(int triidx) {
	if (triangletypes[triidx] == TypeTri::T1)
		triangletypes[triidx] = TypeTri::T1_irregular;
	else if (triangletypes[triidx] == TypeTri::T2)
		triangletypes[triidx] = TypeTri::T2_irregular;
	else if (triangletypes[triidx] == TypeTri::T0) {
		;
	}
}