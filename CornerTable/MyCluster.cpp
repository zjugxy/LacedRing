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

        dualnodes[id1].version = -1;
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

    LacedWireGenerator(mesh);


}

void MyCluster::LacedWireGenerator(const MyMesh& mesh)
{
    vertexsets.resize(dualnodes.size());
    //先生成vertex sets
    int dualid = 0;
    for (auto& elem : dualnodes) {
        //若node 有效
        if (elem.version != -1&&elem.faces.size()!=0) {
            for (auto& ver : elem.vertices) {
                int cnt = VertexOccurCnt(ver, dualid);
                if (cnt >= 2)vertexsets[dualid].boundset.insert(ver);
                if (cnt > 2)vertexsets[dualid].cornerset.insert(ver);
            }
            //被包围了，所以没有corner
            if ((vertexsets[dualid].boundset.size() != 0) && (vertexsets[dualid].cornerset.size()==0)) {
                vertexsets[dualid].surrounded = 1;//被包围
                if (adj_graph[dualid].size() != 1) {
                    std::cout << "error may occur in surrounded mark" << std::endl;
                }
                int adjidx = adj_graph[dualid].begin()->first;
                vertexsets[adjidx].surrounded = 2;//包围了
            }

            for (auto& face : elem.faces)
                vertexsets[dualid].faces.insert(face);
        }
        dualid++;
    }

    lacewires.resize(dualnodes.size());


    for (int i = 0; i < vertexsets.size(); i++) {
        if (vertexsets[i].boundset.size() != 0) {
            LaceWires temp;
            temp.dualid = i;
            {
                if(i==2590)
                    std::cout << "debug　" << i << std::endl;
                if (vertexsets[i].surrounded != 0) {
                    lateaddress.push_back(i);
                    continue;
                }

                LoadVertexset2Wire(temp, vertexsets[i], mesh);
            }
            lacewires[i] = (temp);
        }
    }


    for (auto& id : lateaddress) {
        if (vertexsets[id].surrounded == 1) {
            if (vertexsets[id].cornerset.size() != 0)std::cout << "error lateaddress sur ==1" << std::endl;
            //先处理lacewires[i]
            int adjid = (adj_graph[id].begin())->first;
            LaceWires temp;
            temp.dualid = id;
            LoadSingleLoop(temp, vertexsets[id], mesh);
            lacewires[id] = (temp);
            lacewires[adjid].wires.push_back(temp.wires[0]);
        }
    }

    for (auto& id : lateaddress) {
        if (vertexsets[id].surrounded == 2) {
            //lacewires[id]其中的wires已经被一部分初始化了
            LaceWires& temp = lacewires[id];
            if (temp.wires.size() == 0)std::cout << "error may occur in wire surround others" << std::endl;
            LoadVertexset2WireVersion2(temp, vertexsets[id], mesh);

        }
    }


}

void MyCluster::PackintoLaceWire(std::vector<ExternalWire>& ewires, std::vector<LaceWire_meshlet>& meshlets, std::map<int, int>& dual2idx)
{
    std::map<std::array<int, 4>, int> wire4idx2id;
    //pack这边要调整
    int extcnt = 0;
    //可能可以在这里处理掉许多corner case
    for(const auto&wireloop:lacewires)
        for (const auto& wire : wireloop.wires) {
            std::array<int, 4> temp = packvec2array(wire);
            if (wire4idx2id.find(temp) == wire4idx2id.end()) {
                ExternalWire ewire;
                for (const auto& elem : wire)ewire.wire.push_back(elem);
                ewires.push_back(ewire);
                wire4idx2id[temp] = extcnt++;
            }
        }

    int meshletcnt = 0;
    for(int i=0;i<dualnodes.size();i++)
        if (dualnodes[i].version != -1) {
            LaceWire_meshlet lwmeshlet;
            if (meshletcnt == 72)
                std::cout << "debug meshlet 72" << std::endl;

            for (const auto& wire : lacewires[i].wires) 
                lwmeshlet.externalwireids.push_back(wire4idx2id[packvec2array(wire)]);
            for (const auto& elem : dualnodes[i].faces)
                lwmeshlet.faces.insert(elem);
            for (const auto& elem : adj_graph[i])
                lwmeshlet.adj_nodes.insert(elem.first);
            for (const auto& elem : dualnodes[i].vertices)
                lwmeshlet.vertexs.insert(elem);

            lwmeshlet.facenum = lwmeshlet.faces.size();
            lwmeshlet.vertexnum = lwmeshlet.vertexs.size();
            meshlets.push_back(lwmeshlet);
            dual2idx[i] = meshletcnt++;
        }
    //处理有包围情况的meshlet的external wire
    for (auto& id : lateaddress) {
        if (vertexsets[id].surrounded == 1) {
            if (vertexsets[id].cornerset.size() != 0)std::cout << "error lateaddress sur ==1" << std::endl;
            //先处理lacewires[i]
            int adjid = (adj_graph[id].begin())->first;
            int meshletid = dual2idx[id];
            int adjmeshletid = dual2idx[adjid];
            meshlets[adjmeshletid].externalwireids.push_back(meshlets[meshletid].externalwireids[0]);
        }
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

void MyCluster::LoadVertexset2Wire(LaceWires& wire, const VertexSet& vertexset, const MyMesh& mesh)
{
    int nodeid = wire.dualid;
    int start = *(vertexset.cornerset.cbegin());
    int pnt = start;
    std::set<int> pre = {};
    for (int i = 0; i < vertexset.cornerset.size(); i++) {
        std::vector<int> onewire;
        //int nextpnt = FindNextPnt(pnt, pre, vertexset, mesh);
        //pre.insert(pnt);

        while (true)
        {
            onewire.push_back(pnt);
            //nextpnt找不到的
            if (pre.size() == (vertexset.boundset.size() - 1)) {
                onewire.push_back(start);
                wire.wires.push_back(onewire);
                return;
            }
            if (pnt == 1500)
                std::cout << "debug" << std::endl;
            int nextpnt = FindNextPnt(pnt, pre, vertexset, mesh);
            pre.insert(pnt);
            pnt = nextpnt;
            //pnt为corner
            if (vertexset.cornerset.find(pnt) != vertexset.cornerset.end()) {
                onewire.push_back(pnt);
                wire.wires.push_back(onewire);
                break;
            }
            else {
                continue;
            }
        }

    }

}

void MyCluster::LoadVertexset2WireVersion2(LaceWires& wire, const VertexSet& vertexset, const MyMesh& mesh)
{
    int nodeid = wire.dualid;
    int start = *(vertexset.cornerset.cbegin());
    int pnt = start;
    std::set<int> pre = {};

    for (auto& elem : wire.wires)
        for (auto& v : elem)
            pre.insert(v);


    for (int i = 0; i < vertexset.cornerset.size(); i++) {
        std::vector<int> onewire;
        //int nextpnt = FindNextPnt(pnt, pre, vertexset, mesh);
        //pre.insert(pnt);

        while (true)
        {
            onewire.push_back(pnt);
            //nextpnt找不到的
            if (pre.size() == (vertexset.boundset.size() - 1)) {
                onewire.push_back(start);
                wire.wires.push_back(onewire);
                return;
            }
            if (pnt == 1897)
                std::cout << "debug" << std::endl;
            int nextpnt = FindNextPnt(pnt, pre, vertexset, mesh);
            pre.insert(pnt);
            pnt = nextpnt;
            //pnt为corner
            if (vertexset.cornerset.find(pnt) != vertexset.cornerset.end()) {
                onewire.push_back(pnt);
                wire.wires.push_back(onewire);
                break;
            }
            else {
                continue;
            }
        }

    }

}

void MyCluster::LoadWiresur2(LaceWires& wire, const VertexSet& vertexset, const MyMesh& mesh)
{

    int nodeid = wire.dualid;
    int start = *(vertexset.cornerset.cbegin());
    int pnt = start;
    std::set<int> pre = {};
    for (auto& singleloop : wire.wires)
        for (auto& elem : singleloop)
            pre.insert(elem);

    for (int i = 0; i < vertexset.cornerset.size(); i++) {
        std::vector<int> onewire;
        //int nextpnt = FindNextPnt(pnt, pre, vertexset, mesh);
        //pre.insert(pnt);

        while (true)
        {
            onewire.push_back(pnt);
            //nextpnt找不到的
            if (pre.size() == (vertexset.boundset.size() - 1)) {
                onewire.push_back(start);
                wire.wires.push_back(onewire);
                return;
            }
            int nextpnt = FindNextPnt(pnt, pre, vertexset, mesh);
            pre.insert(pnt);
            pnt = nextpnt;
            //pnt为corner
            if (vertexset.cornerset.find(pnt) != vertexset.cornerset.end()) {
                onewire.push_back(pnt);
                wire.wires.push_back(onewire);
                break;
            }
            else {
                continue;
            }
        }

    }
}

void MyCluster::LoadSingleLoop(LaceWires& wire, const VertexSet& vertexset, const MyMesh& mesh)
{
    int start = *vertexset.boundset.begin();
    int pnt = start;
    std::vector<int> onewire;
    std::set<int> pre = {};
    while (true)
    {
        onewire.push_back(pnt);
        //nextpnt找不到的
        if (pre.size() == (vertexset.boundset.size() - 1)) {
            onewire.push_back(start);
            wire.wires.push_back(onewire);
            return;
        }
        int nextpnt = FindNextPnt(pnt, pre, vertexset, mesh);
        pre.insert(pnt);
        pnt = nextpnt;
        //pnt为corner
    }
    std::cout << "error in single loop" << std::endl;

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
    //return Efit + alphas * Eshape ;
    double addelem = dualnodes[id1].num_face + dualnodes[id2].num_face;
    return Efit + alphas * Eshape  + (addelem/max_face)*(addelem/max_face);
}

int MyCluster::VertexOccurCnt(int vertexid, int dualnodeid) const
{
    int cnt = 1;
    for (auto& neigh : adj_graph[dualnodeid]) 
        if (dualnodes[neigh.first].version != -1) {
            //邻居有效
            if (dualnodes[neigh.first].vertices.find(vertexid) != dualnodes[neigh.first].vertices.end()) {
                cnt++;
            }
        }

    return cnt;
}

bool MyCluster::VerticeInsertCheck(int id1, int id2,std::unordered_set<uint>& test)
{
    for (auto& elem : dualnodes[id1].vertices)test.insert(elem);
    for (auto& elem : dualnodes[id2].vertices)test.insert(elem);
    if (test.size() > max_vertex)return false;
    return true;
}

int MyCluster::FindNextPnt(int pnt, std::set<int>& pre,const VertexSet& vertexset,const MyMesh& mesh) const
{
    auto vhl = mesh.vertex_handle(pnt);
    

    //for (MyMesh::VertexVertexIter vv_it = mesh.cvv_begin(vhl); vv_it != mesh.cvv_end(vhl); ++vv_it) {
    //    int idx = vv_it->idx();
    //    if (  (pre.find(idx) == pre.end()) 
    //        && (vertexset.boundset.find(idx) != vertexset.boundset.end()) ) {
    //        return idx;
    //    }
    //}
    for (auto ve_it = mesh.cve_begin(vhl); ve_it != mesh.cve_end(vhl); ++ve_it) {
        int fidx0 = ve_it->h0().face().idx();
        int fidx1 = ve_it->h1().face().idx();
        
        int count0 = vertexset.faces.count(fidx0);
        int count1 = vertexset.faces.count(fidx1);
        //不是边界边
        if ((count0 == 1) && (count1 == 1))continue;
        else if ((count0 == 0) && (count1 == 0))continue;
        else {
            int temp = ve_it->v0().idx();
            int idx = (temp == pnt) ? ve_it->v1().idx() : temp;
            if (pre.find(idx) == pre.end())return idx;
            else continue;
        }

    }

    std::cout << "error in find next pnt" << std::endl;
    return -1;
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

std::array<int, 4> MyCluster::packvec2array(const std::vector<int>& wire)
{
    std::array<int, 4> temp;
    temp[0] = wire[0]; temp[1] = wire[1];
    int lastid = wire.size() - 1;
    temp[3] = wire[lastid]; temp[2] = wire[lastid - 1];
    return temp;
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