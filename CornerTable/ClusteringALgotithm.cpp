#include "ClusteringALgotithm.h"
#include <queue>
#include <algorithm>
#define MAXFACENUMINMESHLET 60


DualNodeType::DualNodeType(int f, const Eigen::Matrix3d& A, const Eigen::Vector3d& b, double c) : A(A), b(b), c(c) {
    faces = std::unordered_set<int>();
    faces.insert(f);
}

DualNodeType DualNodeType::operator+(const DualNodeType& rhs) const {
    DualNodeType res;
    res.faces = faces;
    for (int face : rhs.faces)
        res.faces.insert(face);
    res.A = A + rhs.A;
    res.b = b + rhs.b;
    res.c = c + rhs.c;
    return res;
}

std::vector<std::vector<int>> ClusteringAlgorithm::Cluster(MyMesh mesh, double th) {
    // V�洢vertices geometry information
    // F�洢Face index information

    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    //igl::readOBJ(input_filename, V, F);

    std::vector<DualNodeType> dual_nodes;
    std::unordered_set<int> valid;
    //raw version
    //for (int i = 0; i < F.rows(); ++i) {
    //    int i0 = F(i, 0), i1 = F(i, 1), i2 = F(i, 2);
    //    DualNodeType dn0(i, V.row(i0).transpose().eval() * V.row(i0), V.row(i0).transpose(), 1);
    //    DualNodeType dn1(i, V.row(i1).transpose().eval() * V.row(i1), V.row(i1).transpose(), 1);
    //    DualNodeType dn2(i, V.row(i2).transpose().eval() * V.row(i2), V.row(i2).transpose(), 1);
    //    dual_nodes.emplace_back(dn0 + dn1 + dn2);
    //    valid.insert(i);
    //}
    //openmesh rewrite version
    for (auto& face : mesh.faces()) {
        std::vector<Eigen::Vector3<double>> vecs;
        int id = face.idx();
        //std::cout << id << std::endl;
        for (const auto& ver : face.vertices()) {
            auto pnt = mesh.point(ver);
            vecs.push_back(Eigen::Vector3<double>{pnt[0], pnt[1], pnt[2]});
        }
        DualNodeType dn0(id, vecs[0] * vecs[0].transpose(), vecs[0], 1);
        DualNodeType dn1(id, vecs[1] * vecs[1].transpose(), vecs[1], 1);
        DualNodeType dn2(id, vecs[2] * vecs[2].transpose(), vecs[2], 1);
        dual_nodes.emplace_back(dn0 + dn1 + dn2);
        valid.insert(id);
    }

    std::vector<std::pair<int, int>> dual_edges;
    GetDualEdges(mesh, dual_edges);

    std::vector<std::vector<int>> dual_adjacency_list(dual_nodes.size());
    //�Ե�һ��Ԫ�ؽ�������
    std::priority_queue<std::pair<double, int>> pq;
    for (int i = 0; i < dual_edges.size(); ++i) {
        const auto& p = dual_edges[i];
        const DualNodeType& n1 = dual_nodes[p.first], & n2 = dual_nodes[p.second];
        double cost = EvaluateCost(n1, n2);
        pq.push({ -cost, i });

        dual_adjacency_list[p.first].push_back(p.second);
        dual_adjacency_list[p.second].push_back(p.first);
    }
    //�ܺϲ���dual edge���ϲ���
    while (!pq.empty()) {
        double cost = -pq.top().first;
        if (cost > th) break;
        auto p = dual_edges[pq.top().second];//pΪҪ�ϲ���dual edge
        pq.pop();
        int nd_id0 = p.first, nd_id1 = p.second;
        if (!valid.count(nd_id0) || !valid.count(nd_id1)) continue;
        if ((dual_nodes[nd_id0].faces.size() + dual_nodes[nd_id1].faces.size()) > MAXFACENUMINMESHLET)continue;
        //�ϲ��������µ�node
        auto new_id = static_cast<int>(dual_nodes.size());
        dual_nodes.emplace_back(dual_nodes[nd_id0] + dual_nodes[nd_id1]);

        std::vector<int> new_neighborhood;
        for (int neighbor : dual_adjacency_list[nd_id0])
            if (valid.count(neighbor) && neighbor != nd_id1) {
                new_neighborhood.push_back(neighbor);
                dual_adjacency_list[neighbor].push_back(new_id);

                double new_cost = EvaluateCost(dual_nodes[nd_id0], dual_nodes[neighbor]);
                pq.push({ -new_cost, static_cast<int>(dual_edges.size()) });
                dual_edges.emplace_back(neighbor, new_id);
            }
        for (int neighbor : dual_adjacency_list[nd_id1])
            if (valid.count(neighbor) && neighbor != nd_id0) {
                new_neighborhood.push_back(neighbor);
                dual_adjacency_list[neighbor].push_back(new_id);

                double new_cost = EvaluateCost(dual_nodes[nd_id1], dual_nodes[neighbor]);
                pq.push({ -new_cost, static_cast<int>(dual_edges.size()) });
                dual_edges.emplace_back(neighbor, new_id);
            }
        dual_adjacency_list.push_back(new_neighborhood);

        valid.insert(new_id);
        valid.erase(nd_id0), valid.erase(nd_id1);
        dual_adjacency_list[nd_id0].clear();
        dual_adjacency_list[nd_id1].clear();
    }

    std::vector<std::vector<int>> res;
    for (auto nd_id : valid) {
        std::vector<int> patch;
        for (auto f : dual_nodes[nd_id].faces) {
            patch.push_back(f);
            //std::cout << f << ' ';
        }
        res.push_back(patch);
        //std::cout << std::endl;
    }

    std::cout << "meshlet generation successfully" << std::endl;
    return res;
}

void ClusteringAlgorithm::GetDualEdges(MyMesh mesh,
    std::vector<std::pair<int, int>>& dual_edges) {
    for (auto& edge : mesh.edges()) {
        auto hf0 = edge.halfedge(0);
        auto hf1 = edge.halfedge(1);
        if (hf0.is_boundary() || hf1.is_boundary())
            continue;
        dual_edges.emplace_back(hf0.face().idx(), hf1.face().idx());
    }
}

double ClusteringAlgorithm::EvaluateCost(const DualNodeType& dn1, const DualNodeType& dn2) {
    Eigen::Matrix3d A = dn1.A + dn2.A;
    Eigen::Vector3d b = dn1.b + dn2.b;
    double c = dn1.c + dn2.c;
    //���Z������ֵ��������
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

    return ((n.transpose().eval() * A).dot(n) + 2 * d * n.transpose().eval().dot(b)) / c + d * d;
}
