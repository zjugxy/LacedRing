#include <iostream>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
typedef OpenMesh::TriMesh_ArrayKernelT<>  MyMesh;
#include <Eigen/Dense>
#include"Viewer/scene.h"
#include"Viewer/viewer.h"
#include"ClusteringALgotithm.h"
#include"MyCluster.h"
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

    //ClusteringAlgorithm clu;
    //auto meshlets = clu.Cluster(mesh, 0.02);

    MyCluster clu(mesh, 64, 126,0.5);
    auto meshlets = clu.oldmeshlets;

    int dd = 0;
    for (const auto& meshlet : meshlets) {
        dd += meshlet.size();
    }
    std::cout << "total face in meshlets " << dd << std::endl;
    std::cout << "total face in mesh" << mesh.n_faces() << std::endl;
    //


    //渲染hamiltoniancycle环，需要设置geometryshader中的faceidtocolor数组大小
    {
        glfwviewer::Viewer myview;
        myview.initGLFW();
        glfwviewer::Scene myscene;

        //myscene.LoadTSMeshlet(mesh, meshlets);
        myscene.LoadCornerTable(mesh);

        /*
        need edit following three line
        */
        myscene.LoadSCMeshlet(mesh, meshlets);
        myview.set(&myscene);
        myview.setMeshshader("SCmeshshader.glsl", "TSfragshader.glsl");

        while (!glfwWindowShouldClose(myview.MYwindow()))
        {
            myview.processinput();
            //myview.RenderuseLR();
            /*
            need edit to ust different meshlets type
            */
            myview.RenderSCML();

            glfwSwapBuffers(myview.MYwindow());
            glfwPollEvents();
        }
        glfwTerminate();
        exit(EXIT_SUCCESS);
    }




    return 0;
}
