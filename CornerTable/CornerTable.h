#pragma once
#include<vector>
#include<array>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include"Viewer/viewermath.h"
#include<list>

typedef OpenMesh::TriMesh_ArrayKernelT<>  MyMesh;

enum TypeTri
{
	T0,
	T1,
	T2,
	T1_irregular,
	T2_irregular
};

class CornerTable
{
public:
	CornerTable(MyMesh);

	//void VTableTest();


	//corner op
	vec3 getVertexValue(int corner)const;
	int Vertex(int corner)const;
	int Tri(int corner)const;
	int Next(int corner)const;
	int Prev(int corner)const;
	int Left(int corner)const;
	int Right(int corner)const;
	int Opp(int corner)const;
	int CornerofVetex(int vertex)const;
	int Swin(int corner)const;

	void TriangleClassification();
	
	//a new version of hamiton cycle with better vertices ring generation
	void BetterHamiltonianCycle(int start);

private:
	bool simplecheck(int triid, std::list<int>::iterator& it)const;
	//
	//ע�͵��ҽ�����ǰ�����α�invadeʱ�������Ӧ��vertex,
	//��ǰ�߰�����
	void HamiltonianCycle(int start);
	int Vertex2Tris(int v0,int v1,int& corner)const;
	void Tri_Ring_Add(int triid);
	void AdjustTritype_AdjecentT0(int triidx);


public:
	/*
	CT basic ��Ա����
	*/
	std::vector<int> Vtable;//Vtable[corner] ���Զ�ȡ��corner��Ӧ��vertex
	std::vector<int> Otable;//���Զ�ȡ��corner��Ӧ��oppo
	std::vector<int> Ctable;//Ctable[vertexindex] = corner
	std::vector<vec3> vertices;

	/*
	temp variables for Ring generate
	*/



	std::vector<int> Mvertex;//mark vertex whether has been visited
	std::vector<int> Mtri;//triangle whether has been visited

	std::vector<int> ringcorner;
	std::vector<int> ringvertex;

	std::vector<int> RingRightCorners;
	std::vector<int> RingLeftCorners;

	std::vector<TypeTri> triangletypes;
	std::vector<int> RenderTriType;
	std::vector<int> T0trilist;

	//���ݽṹ������Ⱦ�����Լ�T0�����֤��rightcorners and leftcorners�ǶԵ�
	std::vector<int> facerightleftwithcolor;
};

