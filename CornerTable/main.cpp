#include <iostream>
 // -------------------- OpenMesh
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include "CornerTable.h"
// ----------------------------------------------------------------------------
#include"Viewer/scene.h"
#include"Viewer/viewer.h"
#include"MyLR.h"
typedef OpenMesh::TriMesh_ArrayKernelT<>  MyMesh;

int main()
{
    MyMesh mesh;

    // generate vertices
    try
    {
        if (!OpenMesh::IO::read_mesh(mesh, "bunnysimple.obj"))
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
    
    //{
    //std::cout << "vertices part" << std::endl;
    //for (const auto& ver : mesh.vertices()) {
    //    std::cout << ver.idx() << ' ';
    //    auto pnt = mesh.point(ver);
    //    std::cout << pnt[0] << ' ' << pnt[1] << ' ' << pnt[2] << std::endl;
    // }
    //}

    //std::cout << "face part" << std::endl;
    //for (const auto& face : mesh.faces()) {
    //    std::cout << face.idx() << ' ';
    //    for (const auto& ver : face.vertices())
    //        std::cout << ver.idx() << ' ';
    //    std::cout << std::endl;
    //}

    

    CornerTable ct(mesh);
    ct.BetterHamiltonianCycle(9);
    ct.TriangleClassification();

    //使用LR对CornerTable进行压缩

    MyLR lr(ct);

    //渲染hamiltoniancycle环，需要设置geometryshader中的faceidtocolor数组大小
    {
        glfwviewer::Viewer myview;
        myview.initGLFW();
        glfwviewer::Scene myscene;
        myscene.LoadCornerTable(ct);
        myscene.LoadLR(lr);

        myview.set(&myscene);
        myview.setshader("vertexshader.glsl", "fragshader.glsl", "geometryshader.glsl");
        myview.setringshader("ringvertex.glsl", "ringfrag.glsl");



        while (!glfwWindowShouldClose(myview.MYwindow()))
        {
            myview.processinput();
            //false时渲染Hamilton
            //true时渲染triangle classification
            // 
             
            myview.RenderuseLR();
            //myview.Renderct(false);
            myview.Renderring();


            glfwSwapBuffers(myview.MYwindow());
            glfwPollEvents();
        }
        glfwTerminate();
        exit(EXIT_SUCCESS);
    }

    return 0;
}

// 运行程序: Ctrl + F5 或调试 >“开始执行(不调试)”菜单
// 调试程序: F5 或调试 >“开始调试”菜单

// 入门使用技巧: 
//   1. 使用解决方案资源管理器窗口添加/管理文件
//   2. 使用团队资源管理器窗口连接到源代码管理
//   3. 使用输出窗口查看生成输出和其他消息
//   4. 使用错误列表窗口查看错误
//   5. 转到“项目”>“添加新项”以创建新的代码文件，或转到“项目”>“添加现有项”以将现有代码文件添加到项目
//   6. 将来，若要再次打开此项目，请转到“文件”>“打开”>“项目”并选择 .sln 文件

