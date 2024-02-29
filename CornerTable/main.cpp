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
#include"Meshlet.h"
bool meshletexist = false;
std::string meshletfilename;
    
#define USINGSC 1



int main()
{
    MyMesh mesh;
    Meshlets meshlets;
    std::string filename = "E:/OpenMesh/Models/horse/horse.obj";
    std::string debugfile = "E:/OpenMesh/Models/horse/horse.meshlets";
    // generate vertices
    meshlets = readVectorFromFile(debugfile);


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



    //meshlets = readVectorFromFile(meshletfilename);
    //if (!meshlets.empty()) {
    //    meshletexist = true;
    //}

    std::cout <<"face num:" << mesh.n_faces() << std::endl<<"vert num:" << mesh.n_vertices() << std::endl <<"halfedge :" << mesh.n_halfedges() << std::endl;
    

        //MyCluster clu(mesh, 64, 126, 0.5);
        //meshlets = clu.oldmeshlets;

        NewCluster nclu(64, 126,mesh,meshlets,debugfile);
        meshlets = nclu.oldmeshlets;


        int flag = 1;
        glfwviewer::Viewer myview;
        myview.initGLFW();
        glfwviewer::Scene myscene;
        myscene.LoadCornerTable(mesh);


        if (flag == 0) {
            myscene.LoadSCMeshlet(mesh, meshlets);
            //myscene.LoadPoints(mesh);
            //myscene.LoadInternalWire(mesh, nlwgen);
            //myscene.LoadLaceWire(mesh, nclu);
            //myscene.LoadNormalLine(nclu);

            //myscene.LoadLines(mesh, meshlets);
            myview.set(&myscene);
            myview.setMeshshader("SCmeshshader.glsl", "TSfragshader.glsl");
            myview.setlineshader("ringvertex.glsl", "ringfrag.glsl");

            while (!glfwWindowShouldClose(myview.MYwindow()))
            {
                myview.processinput();
                myview.RenderSCML();
                //myview.RenderPoint();
                //myview.RenderNormalLine();
                //myview.RenderWireLine();
                //myview.RenderNormalLine();
                //myview.RenderInterWire();
                glfwSwapBuffers(myview.MYwindow());
                glfwPollEvents();
            }
            glfwTerminate();
            exit(EXIT_SUCCESS);
        }
        else if(flag ==1){
            NewLWGenerator nlwgen(nclu);
            nlwgen.FUNC(mesh);
            myscene.LoadSimpleWireMeshlet(nlwgen);
            //myscene.LoadInternalWire(mesh, nlwgen);
            myscene.LoadLaceWire(mesh, nclu);
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
        else if(flag == 2){
            NewLWGenerator nlwgen(nclu);
            nlwgen.FUNC(mesh);
            myscene.LoadGPULW(nlwgen);
            myscene.LoadInternalWire(mesh, nlwgen);
            myscene.LoadLaceWire(mesh, nclu);
            myview.set(&myscene);
            myview.setMeshshader("Lacemeshshader.glsl", "TSfragshader.glsl");
            myview.setlineshader("ringvertex.glsl", "ringfrag.glsl");

            while (!glfwWindowShouldClose(myview.MYwindow()))
            {
                myview.processinput();
                myview.RenderFINALLWML();
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
