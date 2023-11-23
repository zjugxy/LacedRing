#include <iostream>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
typedef OpenMesh::TriMesh_ArrayKernelT<>  MyMesh;
#include <Eigen/Dense>
#include"Viewer/scene.h"
#include"Viewer/viewer.h"
#include"ClusteringALgotithm.h"
#include"MyCluster.h"
#include"LaceWireGenerator.h"
bool meshletexist = false;
std::string meshletfilename;

int main()
{
    MyMesh mesh;
    Meshlets meshlets;
    std::string filename = "E:/OpenMesh/sphere.obj";
    // generate vertices
    try
    {
        if (!OpenMesh::IO::read_mesh(mesh, filename))
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

    meshletfilename = changeFileExtension(filename, "meshlet");

    //meshlets = readVectorFromFile(meshletfilename);
    //if (!meshlets.empty()) {
    //    meshletexist = true;
    //}

    std::cout <<"face num:" << mesh.n_faces() << std::endl<<"vert num:" << mesh.n_vertices() << std::endl <<"halfedge :" << mesh.n_halfedges() << std::endl;
    
   // if (meshletexist == false) {
        MyCluster clu(mesh, 64, 126, 0.5);
        meshlets = clu.oldmeshlets;
        writeVectorToFile(meshlets, meshletfilename);

        LaceWireGenerator lwgen;
        clu.PackintoLaceWire(lwgen.Ewires, lwgen.meshlets, lwgen.Dual2idx);

        lwgen.InternalWireGeneraotr(mesh);

   // }

    //int dd = 0;
    //for (const auto& meshlet : meshlets) {
    //    dd += meshlet.size();
    //}
    //std::cout << "total face in meshlets " << dd << std::endl;
    //std::cout << "total face in mesh" << mesh.n_faces() << std::endl;
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
        myscene.LoadLaceWire(mesh,clu);
        //myscene.LoadLacePoint(mesh, clu);
        myscene.LoadSCMeshlet(mesh, meshlets);

        myscene.LoadSimpleWireMeshlet(lwgen);

        myscene.LoadInternalWire(mesh, lwgen);


        myview.set(&myscene);
        myview.setMeshshader("SimpleLaceWiremeshshader.glsl", "TSfragshader.glsl");
        myview.setlineshader("ringvertex.glsl", "ringfrag.glsl");
        while (!glfwWindowShouldClose(myview.MYwindow()))
        {
            myview.processinput();
            //myview.RenderuseLR();
            /*
            need edit to ust different meshlets type
            */
            myview.RenderSCML();
            //myview.RenderWirePnt();
            //myview.RenderSWML();

            myview.RenderWireLine();
            myview.RenderInterWire();


            glfwSwapBuffers(myview.MYwindow());
            glfwPollEvents();
        }
        glfwTerminate();
        exit(EXIT_SUCCESS);
    }

    return 0;
}
