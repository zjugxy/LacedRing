#pragma once
#include<../NewLWGenerator.h>

class VertexEncoder
{
private:
	float xerror, yerror, zerror;

public:
	//����packinter packexter��Ҫ�Ǽ�������xnum,ynum,znum
	//�� geobitflowpack�� datatocompress -> ѹ��
	void CalculateXYZnum(PackGEO& tempgeo);

	void CalculateCornerGeobit(
		const std::vector<uint>& corneridxs,
		const MyMesh& mesh,
		Eigen::MatrixXd& dequanmatrix,
		Eigen::Vector3d& tranvec
	);



};

