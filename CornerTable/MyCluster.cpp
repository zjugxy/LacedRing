#include "MyCluster.h"
#define REVERSEPI 0.31830988618

void CalculateInternalCircle(Eigen::Vector3d& v1, Eigen::Vector3d& v2, Eigen::Vector3d& v3, double& wtri, double& diameter);

MyCluster::MyCluster(const MyMesh& mesh, int maxvertex, int maxface,double alpha_shape=1.0)
{
    max_vertex = maxvertex;
    max_face = maxface;
    alphas = alpha_shape;
    //初始化 adjgraph dualnodes pq
    BuildDualNodes(mesh);

    while (!pq.empty()) {
        //check是否有效
        double cost = pq.top().first;
        auto checkelem = pq.top().second;
        pq.pop();
        int id0 = checkelem[0], version0 = checkelem[1];
        int id1 = checkelem[2], version1 = checkelem[3];

        if ((version0 != dualnodes[id0].version) || (version1 != dualnodes[id1].version))continue;
        if (dualnodes[id0].num_face + dualnodes[id1].num_face > max_face)continue;
        //check vertices
        std::unordered_set<uint> testset;
        if (VerticeInsertCheck(id0, id1,testset) == false)continue;
        //

        // id1 --> id0  dualnode[id1]合并进入dualnode[id0]
        //version == -1 表示该dualnode已经被吞并

        dualnodes[id0].version += 1;
        //dualnodes[id0].costelem = dualnodes[id0].costelem + dualnodes[id1].costelem;
        dualnodes[id0].costelem = RefreshElem(dualnodes[id0].costelem, dualnodes[id1].costelem,id0,id1);

        dualnodes[id0].num_face += dualnodes[id1].num_face;
        for (auto& elem : dualnodes[id1].faces)dualnodes[id0].faces.push_back(elem);
        dualnodes[id0].vertices = testset;
        dualnodes[id0].num_vertex = testset.size();

        double length01 = adj_graph[id0][id1];

        for(auto&elem:adj_graph[id1])
            if (dualnodes[elem.first].version != -1 && elem.first != id0) {

                //adj_graph[elem].insert(id0);
                //adj_graph[elem].erase(id1);
                //adj_graph[id0].insert(elem);
                //id1 id0 同时与 elem.first相连时
                if (adj_graph[id0].find(elem.first) != adj_graph[id0].end()) {
                    adj_graph[id0][elem.first] += length01;
                    adj_graph[elem.first][id0] = adj_graph[id0][elem.first];
                    adj_graph[elem.first].erase(id1);
                }
                else {
                    adj_graph[id0][elem.first] = length01;
                    adj_graph[elem.first][id0] = length01;
                    adj_graph[elem.first].erase(id1);
                }
            }
        adj_graph[id0].erase(id1);

        dualnodes[id1].version == -1;
        dualnodes[id1].vertices.clear(); 
        dualnodes[id1].faces.clear();
        dualnodes[id1].num_face = 0;
        dualnodes[id1].num_vertex = 0;
        dualnodes[id1].costelem = EvaluateElem();

        adj_graph[id1].clear();

        for (auto& elem : adj_graph[id0])
            if ((dualnodes[elem.first].version != -1) &&
                ((dualnodes[elem.first].num_face + dualnodes[id0].num_face) <= max_face)) {
                double tempcost = EvaluateCost(elem.first, id0);
                pq.push(std::pair<double, std::array<int, 4>>{
                    -1*tempcost, { id0,dualnodes[id0].version,elem.first,dualnodes[elem.first].version }
                });
            }
    }

    for (auto& elem : dualnodes) {
        if (elem.version != -1)
            oldmeshlets.push_back(elem.faces);
    }


}

void MyCluster::BuildDualNodes(const MyMesh& mesh)
{
    dualnodes.resize(mesh.n_faces());
    adj_graph.resize(mesh.n_faces());


    for (const auto& face : mesh.faces()) {
        std::vector<Eigen::Vector3<double>> vecs;
        int id = face.idx();
        MyDualNode node;
        for (const auto& ver : face.vertices()) {
            node.vertices.insert(ver.idx());
            auto pnt = mesh.point(ver);
            vecs.push_back(Eigen::Vector3<double>{pnt[0], pnt[1], pnt[2]});
        }
        EvaluateElem dn0(vecs[0] * vecs[0].transpose(), vecs[0], 1);
        EvaluateElem dn1(vecs[1] * vecs[1].transpose(), vecs[1], 1);
        EvaluateElem dn2(vecs[2] * vecs[2].transpose(), vecs[2], 1);
        node.costelem = dn0 + dn1 + dn2;

        CalculateInternalCircle(vecs[0], vecs[1], vecs[2], node.costelem.area, node.costelem.diameter);
        node.costelem.gamma = (node.costelem.diameter * node.costelem.diameter) * (REVERSEPI * 0.25)/node.costelem.area;
        node.version = 0;
        node.num_face = 1;
        node.num_vertex = 3;
        node.faces.push_back(id);
        dualnodes[id] = node;

        for (const auto& ff_iter : face.faces())
            if (!ff_iter.is_boundary()) {
                //计算两个
                adj_graph[id].insert(std::pair<int, double>{ff_iter.idx(), 0.0f});
            }
        
    }

    for (const auto& edge : mesh.edges()) {
        if (edge.is_boundary())continue;
        //此时id0,id1都是faceid 和 dualnodeid
        int id0 = edge.halfedge(0).face().idx();
        int id1 = edge.halfedge(1).face().idx();

        auto ver0 = mesh.point(edge.v0());
        auto ver1 = mesh.point(edge.v1());
        double length = (ver0 - ver1).length();

        adj_graph[id0][id1] = length;
        adj_graph[id1][id0] = length;

        double cost = EvaluateCost(id0, id1);
        pq.push(std::pair<double, std::array<int, 4>>{
            -1*cost, { id0,dualnodes[id0].version,id1,dualnodes[id1].version }
        });
    }



}

double MyCluster::EvaluateCost(int id1, int id2)
{
    auto dn1 = dualnodes[id1].costelem;
    auto dn2 = dualnodes[id2].costelem;
    Eigen::Matrix3d A = dn1.A + dn2.A;
    Eigen::Vector3d b = dn1.b + dn2.b;
    double c = dn1.c + dn2.c;
    //Eshape
    double newdiameter = dn1.diameter + dn2.diameter - 2 * adj_graph[id1][id2];
    double newgamma = newdiameter * newdiameter * (0.25 * REVERSEPI) / (dn1.area + dn2.area);
    double Eshape = (newgamma - std::max(dn1.gamma, dn2.gamma)) / newgamma;

    //求解Z的特征值特征向量
    Eigen::Matrix3d Z = A - b * b.transpose().eval() / c;
    Eigen::EigenSolver<Eigen::Matrix3d> eigen_solver(Z);
    std::vector<double> eigenvalues(3);
    std::vector<Eigen::Vector3d> eigenvectors(3);

    for (int i = 0; i < 3; i++) {
        eigenvalues[i] = eigen_solver.eigenvalues()[i].real();
        Eigen::Vector3cd vcd = eigen_solver.eigenvectors().col(i);
        for (int j = 0; j < 3; ++j)
            eigenvectors[i](j) = vcd(j).real();
        eigenvectors[i].normalize();
    }

    for (int i = 0; i < 3; ++i)
        for (int j = i + 1; j < 3; ++j)
            if (eigenvalues[i] > eigenvalues[j]) {
                std::swap(eigenvalues[i], eigenvalues[j]);
                std::swap(eigenvectors[i], eigenvectors[j]);
            }

    Eigen::Vector3d n = eigenvectors[0];
    double d = -n.transpose().dot(b) / c;

    double Efit =  ((n.transpose().eval() * A).dot(n) + 2 * d * n.transpose().eval().dot(b)) / c + d * d ;
    return Efit + alphas * Eshape;
}

bool MyCluster::VerticeInsertCheck(int id1, int id2,std::unordered_set<uint>& test)
{
    for (auto& elem : dualnodes[id1].vertices)test.insert(elem);
    for (auto& elem : dualnodes[id2].vertices)test.insert(elem);
    if (test.size() > max_vertex)return false;
    return true;
}

EvaluateElem EvaluateElem::operator+(const EvaluateElem& rhs) const
{
    EvaluateElem res;
    res.A = A + rhs.A;
    res.b = b + rhs.b;
    res.c = c + rhs.c;
    return res;
}

EvaluateElem MyCluster::RefreshElem(const EvaluateElem& elem1, const EvaluateElem& elem2, int id1, int id2)
{
    EvaluateElem res;
    res.A = elem1.A + elem2.A;
    res.b = elem1.b + elem2.b;
    res.c = elem1.c + elem2.c;

    res.area = elem1.area + elem2.area;
    res.diameter = elem1.diameter + elem2.diameter - 2 * adj_graph[id1][id2];
    res.gamma = res.diameter * res.diameter * (0.25 * REVERSEPI) / (res.area);

    return res;
}


void CalculateInternalCircle(Eigen::Vector3d& v1, Eigen::Vector3d& v2, Eigen::Vector3d& v3, double& wtri, double& diameter) {
    Eigen::Vector3d edge1 = v2 - v1;
    Eigen::Vector3d edge2 = v3 - v1;
    Eigen::Vector3d edge3 = v2 - v3;
    // 计算法向量
    Eigen::Vector3d normal = edge1.cross(edge2);

    // 计算外接圆半径
    wtri = normal.norm() * 0.5;
    double circumcircleRadius = edge1.norm() * edge2.norm() * edge3.norm() / (4 * wtri);
    diameter = 2 * circumcircleRadius;
    return;
}