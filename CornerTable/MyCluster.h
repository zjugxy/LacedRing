#pragma once
#include <unordered_set>
#include<unordered_map>
#include <Eigen/Dense>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include<queue>
typedef OpenMesh::TriMesh_ArrayKernelT<>  MyMesh;

struct EvaluateElem
{
	//最小二乘法项
	Eigen::Matrix3d A;
	Eigen::Vector3d b;
	double c=0.0;

	double area = 0.0;
	double diameter = 0.0;
	double gamma = 0.0;

	EvaluateElem() {};
	EvaluateElem(const Eigen::Matrix3d& A, const Eigen::Vector3d& b, double c) : A(A), b(b), c(c) {}
	EvaluateElem operator+(const EvaluateElem& rhs)const;

};


struct MyDualNode {
	EvaluateElem costelem;
	std::unordered_set<uint> vertices;
	std::vector<int> faces;
	int version = 0;
	int num_face;
	int num_vertex;
};


class MyCluster
{
public:
	std::vector<std::unordered_map<int,double>> adj_graph;//dualnodeid --> dualnodeadj id
	std::vector<MyDualNode> dualnodes;
	std::vector<std::vector<int>> oldmeshlets;
	//double --> cost   cost越小，越排在前面
	//int , 4     --> node0 node0version node1 node1version
	std::priority_queue < std::pair<double, std::array<int, 4>>> pq;

	int max_vertex;
	int max_face;
	double alphas;

public:
	MyCluster(const MyMesh& mesh, int maxvertex, int maxface,double alpha_shape);

private:
	void BuildDualNodes(const MyMesh& mesh);

	double EvaluateCost(int id1, int id2);

	bool VerticeInsertCheck(int id1, int id2,std::unordered_set<uint>& testset);

	EvaluateElem RefreshElem(const EvaluateElem& elem1, const EvaluateElem& elem2
	,int id1,int id2);

};

