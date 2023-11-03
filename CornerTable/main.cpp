#include <iostream>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
typedef OpenMesh::TriMesh_ArrayKernelT<>  MyMesh;
#include <Eigen/Dense>
#include"Viewer/scene.h"
#include"Viewer/viewer.h"
#include"ClusteringALgotithm.h"
int main()
{
    MyMesh mesh;

    // generate vertices
    try
    {
        if (!OpenMesh::IO::read_mesh(mesh, "E:/OpenMesh/sphere.obj"))
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


    //渲染hamiltoniancycle环，需要设置geometryshader中的faceidtocolor数组大小
    {
        glfwviewer::Viewer myview;
        myview.initGLFW();
        glfwviewer::Scene myscene;
        myscene.LoadMeshlet(mesh, meshlets);

        myview.set(&myscene);
        myview.setshader("vertexshader.glsl", "fragshader.glsl", "geometryshader.glsl");

        while (!glfwWindowShouldClose(myview.MYwindow()))
        {
            myview.processinput();
            //myview.RenderuseLR();
            myview.RenderML();

            glfwSwapBuffers(myview.MYwindow());
            glfwPollEvents();
        }
        glfwTerminate();
        exit(EXIT_SUCCESS);
    }




    return 0;
}
