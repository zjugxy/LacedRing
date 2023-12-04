﻿#include <iostream>
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

#define USINGSC 1

int main()
{
    MyMesh mesh;
    Meshlets meshlets;
    std::string filename = "E:/OpenMesh/horse.obj";
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
    

        MyCluster clu(mesh, 64, 126, 0.5);
        meshlets = clu.oldmeshlets;
        LaceWireGenerator lwgen;
        clu.PackintoLaceWire(lwgen.Ewires, lwgen.meshlets, lwgen.Dual2idx);
        lwgen.InternalWireGeneraotr(mesh);

    
        glfwviewer::Viewer myview;
        myview.initGLFW();
        glfwviewer::Scene myscene;
        myscene.LoadCornerTable(mesh);

        if (0) {
            myscene.LoadSCMeshlet(mesh, meshlets);
            //myscene.LoadInternalWire(mesh, lwgen);
            //myscene.LoadLaceWire(mesh, clu);
            myscene.LoadLaceWire(lwgen, mesh);
            myview.set(&myscene);
            myview.setMeshshader("SCmeshshader.glsl", "TSfragshader.glsl");
            myview.setlineshader("ringvertex.glsl", "ringfrag.glsl");

            while (!glfwWindowShouldClose(myview.MYwindow()))
            {
                myview.processinput();
                myview.RenderSCML();
                myview.RenderWireLine();
                //myview.RenderInterWire();
                glfwSwapBuffers(myview.MYwindow());
                glfwPollEvents();
            }
            glfwTerminate();
            exit(EXIT_SUCCESS);
        }
        else {

            myscene.LoadSimpleWireMeshlet(lwgen);
            //myscene.LoadInternalWire(mesh, lwgen);
            //myscene.LoadLaceWire(mesh, clu);
            myscene.LoadLaceWire(lwgen, mesh);
            myview.set(&myscene);
            myview.setMeshshader("SimpleLaceWiremeshshader.glsl", "TSfragshader.glsl");
            myview.setlineshader("ringvertex.glsl", "ringfrag.glsl");

            while (!glfwWindowShouldClose(myview.MYwindow()))
            {
                myview.processinput();
                myview.RenderSWML();
                myview.RenderWireLine();
                //myview.RenderInterWire();
                glfwSwapBuffers(myview.MYwindow());
                glfwPollEvents();
            }
            glfwTerminate();
            exit(EXIT_SUCCESS);
        }

    

    return 0;
}
