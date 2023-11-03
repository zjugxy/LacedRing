#include <iostream>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
typedef OpenMesh::TriMesh_ArrayKernelT<>  MyMesh;
#include <Eigen/Dense>

#include"ClusteringALgotithm.h"
int main()
{
    MyMesh mesh;

    // generate vertices
    try
    {
        if (!OpenMesh::IO::read_mesh(mesh, "E:/OpenMesh/wtbunny.obj"))
        {
            std::cerr << "Cannot write mesh to file 'output.off'" << std::endl;
            return 1;
        }
    }
    catch (std::exception& x)
    {
        std::cerr << x.what() << std::endl;
        return 1;
    }

    std::cout <<"face num:" << mesh.n_faces() << std::endl<<"vert num:" << mesh.n_vertices() << std::endl <<"halfedge :" << mesh.n_halfedges() << std::endl;
    
    using Eigen::MatrixXd;

    ClusteringAlgorithm clu;
    auto meshlets = clu.Cluster(mesh, 0.02);

    //通过渲染验证对不对




    return 0;
}
