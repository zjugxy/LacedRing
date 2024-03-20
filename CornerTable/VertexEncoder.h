#pragma once
#include<../NewLWGenerator.h>

class VertexEncoder
{
private:
	float xerror, yerror, zerror;

public:
	//对于packinter packexter主要是计算它的xnum,ynum,znum
	//在 geobitflowpack中 datatocompress -> 压缩
	void CalculateXYZnum(PackGEO& tempgeo);

	void CalculateCornerGeobit(
		const std::vector<uint>& corneridxs,
		const MyMesh& mesh,
		Eigen::MatrixXd& dequanmatrix,
		Eigen::Vector3d& tranvec
	);



};

