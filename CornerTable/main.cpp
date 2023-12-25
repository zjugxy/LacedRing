#include <iostream>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
typedef OpenMesh::TriMesh_ArrayKernelT<>  MyMesh;
#include <Eigen/Dense>
#include"Viewer/scene.h"
#include"Viewer/viewer.h"
#include"ClusteringALgotithm.h"
#include"NewCluster.h"
#include"NewLWGenerator.h"
#include"MeshFile.h"

#define USINGSC 1

int main()
{
    MyMesh mesh;

    std::string filename = "E:/OpenMesh/Models/armadillo/armadillo.obj";
    // generate vertices

    MeshFile fileaddressor;
    bool fileexist = fileaddressor.ImportFile(filename);

    if (fileexist == false) {
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

        std::cout << "face num:" << mesh.n_faces() << std::endl << "vert num:" << mesh.n_vertices() << std::endl << "halfedge :" << mesh.n_halfedges() << std::endl;


        NewCluster nclu(64, 126, mesh);

        NewLWGenerator nlwgen(nclu);
        nlwgen.FUNC(mesh);
        nlwgen.ExportFile(filename);

        glfwviewer::Scene myscene;
        myscene.LoadMesh(mesh,filename);
        


        std::cout << "file contrust down" << std::endl;
        std::cout << "reopen program" << std::endl;
        return 0;
    }



        glfwviewer::Viewer myview;
        myview.initGLFW();
        glfwviewer::Scene myscene;
        //myscene.LoadCornerTable(mesh);

        //myscene.LoadGPULW(nlwgen);
        //myscene.LoadInternalWire(mesh, nlwgen);
        //myscene.LoadLaceWire(mesh, nclu);
        myscene.LoadGPULW(fileaddressor);
        myview.set(&myscene);
        myview.setMeshshader("Lacemeshshader.glsl", "TSfragshader.glsl");
        myview.setlineshader("ringvertex.glsl", "ringfrag.glsl");
        
        while (!glfwWindowShouldClose(myview.MYwindow()))
        {
            myview.processinput();
            myview.RenderFINALLWML();
            //myview.RenderWireLine();
            //myview.RenderInterWire();
            glfwSwapBuffers(myview.MYwindow());
            glfwPollEvents();
        }
        glfwTerminate();
        exit(EXIT_SUCCESS);

    

    return 0;
}
