#pragma once
#include<vector>
#include<unordered_set>
#include"NewCluster.h"



using uchar = unsigned char;
using uint = unsigned int;

struct MeshletBody {
	std::vector<uchar>left;
	std::vector<uchar> right;
	std::vector<uchar> irregular;
	std::vector<uint> vertex;
};

struct LWMeshlet
{
	MeshletBody InternalLW;
	std::vector<uint> ewireidx;
	std::vector<bool> reverse;
};

struct LWbuilder
{
	std::unordered_set<uint> vertices;
	std::unordered_set<uint> faces;

	std::unordered_set<uint> cornerset;
	std::unordered_set<uint> borderset;
};

struct MeshletDes
{
	uchar ewirenum;
	uchar color[3];

	uchar irrnum;
	uchar numvertex;
	uchar numinver;
	uchar ingeostart;

	std::vector<uchar> numexver;

	uint ingeolocation;// 6+ 26
	uint inconlocation;
	std::vector<uint> exgeolocation;//1+5+26
	std::vector<uint> exconlocation;
};


struct Box2d {
	float miny = 0;
	float minz = 0;
	float maxy = 0;
	float maxz = 0;

	void refresh(float yvalue, float zvalue);
	float getSumLength();
};

struct UnitBox {
	
	float minx = 0;
	float miny = 0;
	float minz = 0;
	float maxx = 0;
	float maxy = 0;
	float maxz = 0;
	void refresh(float xvalue,float yvalue, float zvalue);

};


struct PackGEO {
	bool needtransform;
	bool isXPlatForm = false;
	std::vector<vec3> nopackgeo;


	bool NeedExtraPnt = false;
	std::vector<uint> extraVertex;//当lace wire的个数 >0 但是又<3时，需要添加相邻的点
	//重新来构建



	Eigen::Vector3d euler;
	Eigen::Vector3d translation;
	Eigen::Vector3d scaleelem;

	uchar rotatx, rotaty, rotatz;
	uchar translatex, translatey, translatez;
	uchar scalex, scaley, scalez;

	uchar xnum, ynum, znum;
	std::vector<Eigen::Vector3d> dataunit;
	std::vector<Eigen::Vector3d> datatocompress;

	std::vector<uint> geobitflow;
};


class NewLWGenerator
{
public:
	std::vector<MeshletBody> gloEwires;
	std::vector<LWbuilder> internalbuilders;
	std::vector<LWMeshlet> targets;
	uint nummeshlet;
	std::vector<std::unordered_map<uint, uchar>> vidxmaps;

	//simple lw
	std::vector<Simple_meshlet> simplemeshlets;
	std::vector<vec4> geoinfo;
	std::vector<uchar> priminfo;
	// gpu lw
	std::vector<uint> DesLoc;
	std::vector<uint> Desinfo;

	std::vector<uint> newintercon;
	std::vector<uint> newextercon;
	std::vector<uchar> intercon;//inter: left,right,irregular
	std::vector<uchar> extercon;//exter: left,right
	std::vector<float> intergeo;
	std::vector<float> extergeo;

	std::vector<PackGEO> packexter;
	std::vector<PackGEO> packinter;
	UnitBox MeshTransBox;
	UnitBox MeshScaleBox;

	std::vector<std::vector<int>> bitnums;

public:
	NewLWGenerator(const NewCluster& nclu);

	void FUNC(const MyMesh& mesh);
	void ExportFile(const std::string& filename);

private:
	void InternalWireGenerator(const MyMesh& mesh);
	void LeftRightGenertor(const MyMesh& mesh);


private:
	void LevelPntGen(const std::unordered_set<uint>& bordset,const std::unordered_set<uint>& verset,
		std::vector<std::unordered_set<uint>>& lpnts,const MyMesh&mesh);
	void PushAtBack(std::vector<std::unordered_set<uint>>& lpnts, const MyMesh& mesh, std::deque<uint>& singlewire, std::unordered_set<uint>& pntstoremove);
	void PushAtFront(std::vector<std::unordered_set<uint>>& lpnts, const MyMesh& mesh, std::deque<uint>& singlewire, std::unordered_set<uint>& pntstoremove);
	void PackSimpleLaceWire(const MyMesh& mesh);



	void fitPlane(const MyMesh&mesh,const std::vector<uint>&verset, Eigen::Vector3d& centroid, Eigen::Vector3d& normal);
	Box2d getYZAxis(const MyMesh& mesh, const std::vector<uint>& verset, Eigen::Vector3d xaxis, Eigen::Vector3d& yaxis, Eigen::Vector3d& zaxis, float anglediff,Eigen::Vector3d centroid);
	void getScaleBox(const MyMesh& mesh, const std::vector<uint>& verset, Eigen::Vector3d centroid, UnitBox& scalebox, Eigen::Vector3d xaxis, Eigen::Vector3d& yaxis, Eigen::Vector3d& zaxis);
	
	void NewVertexQuantization(const MyMesh& mesh);

	void PackNoPlane(PackGEO& tempgeo, const MyMesh& mesh, const MeshletBody& ewire);

	void PackGPULW(const MyMesh& mesh);

	void BitSortGen(int highestnum,int minvalue);
	void VertexBitGen(const MyMesh& mesh,float errorpercent);

	std::array<uchar, 4> DiscomposeUint(uint value);
	float PackChar4Float(uchar c0, uchar c1, uchar c2, uchar c3);
	uint PackChar4Uint(uchar c0, uchar c1, uchar c2, uchar c3);
	void getStartAxis(Eigen::Vector3d& yaxis, Eigen::Vector3d& zaxis, const Eigen::Vector3d& normal);

	void NewSimpleCheck(PackGEO& tempgeo, const MyMesh& mesh, const std::vector<uint>& verset);

	void CheckByDequantize(const MyMesh& mesh);

	void FinalScaleGen(const MyMesh& mesh);

	void UcharRadGen();
	void TranslatePack(const MyMesh& mesh, UnitBox TranslateBox);


	uchar myclamp(float value, uchar lowbound, uchar highbound);
	void LimitEigen(Eigen::Vector3d& vec);


	Eigen::Vector3d cutfloatByBit(Eigen::Vector3d rawvalue, uint xnum, uint ynum, uint znum);

	void GEObitflowPack();

	void InsertGeoValue(std::vector<uint>& bitflow, uint startidx, uint length, float rawvalue);

	void printBits(uint32_t value);


	Eigen::Vector3d ReadData(const std::vector<uint>& geobits, uint startidx, uint xnum, uint ynum, uint znum);
	float ReadFloat(const std::vector<uint>& geobits, uint startidx, uint num);

	void PackExtraPnts(PackGEO& tempgeo, const MyMesh& mesh,const std::vector<uint>& vertices);


	void ParseTempGeo(const PackGEO& tempgeo, Eigen::MatrixXd& dequanmatrix, Eigen::Vector3d& tranvec);
};

