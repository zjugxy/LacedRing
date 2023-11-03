//
// Created by Haowen Xu on 2018/6/6.
//

#ifndef HIERARCHICALFACECLUSTERING_CLUSTERINGALGORITHM_H
#define HIERARCHICALFACECLUSTERING_CLUSTERINGALGORITHM_H

#include <unordered_set>
#include <Eigen/Dense>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
typedef OpenMesh::TriMesh_ArrayKernelT<>  MyMesh;

struct DualNodeType {
	std::unordered_set<int> faces;
	Eigen::Matrix3d A;
	Eigen::Vector3d b;
	double c;

	DualNodeType() = default;
	DualNodeType(int f, const Eigen::Matrix3d& A, const Eigen::Vector3d& b, double c);

	DualNodeType operator+(const DualNodeType& rhs) const;
};

class ClusteringAlgorithm {
public:
	// constructors
	ClusteringAlgorithm() = default;

	std::vector<std::vector<int>> Cluster(MyMesh mesh, double th);

private:

	void GetDualEdges(MyMesh mesh, std::vector<std::pair<int, int>>& dual_edges);
	double EvaluateCost(const DualNodeType& dn1, const DualNodeType& dn2);
};

#endif //HIERARCHICALFACECLUSTERING_CLUSTERINGALGORITHM_H