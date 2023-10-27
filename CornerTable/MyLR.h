#pragma once
#include"CornerTable.h"

#define NOPTR 0xfffffffe
#define ISPTR 0xffffffff

#define SETNORMAL 0x0
#define SETIRREGULAR 0x1

#define ISIRREGULAR 0x00000001

struct VOleft
{
	int vleft;
	int v0opp;
	int v2opp;
};
struct VOright
{
	int vright;
	int v4opp;
	int v6opp;
};


class MyLR {
public:
	std::vector<vec3> vertices;//vertices的排列顺序取决于corner vertex
	//right left的最末尾1bit
	std::vector<int> RightTable;
	std::vector<int> LeftTable;

	std::vector<int> T0Vtable;
	std::vector<int> T0Otable;
	std::vector<int> Ctable;//???

	std::vector<VOleft> VOLlist;
	std::vector<VOright> VORlist;

	std::vector<int> facetypeid;

	std::vector<int> EBO_check;
	bool first_generate = true;
	int mring, mi, m;//mregular miregular mtotal


public:

	MyLR(const CornerTable& ct);//该构造函数即为压缩方法

	void RenderByLRCHECK();//直接使用LR 进行渲染
	void GenerateEBO();
	void Decompress();//解压缩为corner table?

	void Addresstri(int idx,TypeTri leftt,TypeTri rightt);
	bool ISirregular(int idx);
	bool IsNormaltri(int idx);

	//if irregular return 1 else normal 0
	bool ParseLeft(int leftcode, int& leftvertexidx);
	bool ParseRight(int rightcode, int& rightvertexidx);
	void checkLR();


	int CornertoVertexidx(int corner);
	vec3 CornertoVertexvalue(int corner);
};