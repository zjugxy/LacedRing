#include<memory>
#include "scene.h"
#include<random>
#define EMPTYWIRE 0xFFFFFFFF
#define WIREEND 0xFF

namespace glfwviewer {
    Scene::Scene():ctptr(nullptr)
    {
    }
    Scene::~Scene()
    {

    }
    void Scene::initialize()
    {
        bdSphere()->center(vec3(0, 0, 0));
        bdSphere()->radius(1);
        bdBox()->mymin(vec3(-1, -1, -1));
        bdBox()->mymax(vec3(1, 1, 1));
    }



    void Scene::LoadCornerTable(const CornerTable& ct)
    {
        //更新mmbox mmsphere
        BoundingBox box;

        if (ct.vertices.empty()) {
            std::cout << "empty ct";
            return;
        }

        box.mymin(glm::vec3(ct.vertices[0][0], ct.vertices[0][1], ct.vertices[0][2]));
        box.mymax(glm::vec3(ct.vertices[0][0], ct.vertices[0][1], ct.vertices[0][2]));

        for (const auto& v : ct.vertices) {
            glm::vec3 p(v[0], v[1], v[2]);
            box.mymin(glm::min(box.mymin(), p));
            box.mymax(glm::max(box.mymax(), p));
        }
        m_bbox = box;
        vec3 center = (box.mymax() + box.mymin()) * 0.5f;
        float radius = glm::distance(box.mymax(), box.mymin())*0.5f;
        m_bsphere.center(center);
        m_bsphere.radius(radius);


        //
        ctptr = &ct;
        SceneRenderCT temp;
        //load scenerenderct

        glGenVertexArrays(1, &temp.VAO);
        //glGenBuffers(1, &ctobj.facecolorVBO);
        glGenBuffers(1, &temp.vertexVBO);
        glGenBuffers(1, &temp.EBO);
        glGenBuffers(1, &temp.facecolorVBO);
        glGenBuffers(1, &temp.lineEBO);



        temp.needrender = true;
        temp.ringneedrender = true;
        ctobj = temp;

        glActiveTexture(GL_TEXTURE0);
        glGenTextures(1, &ctobj.texture);
        glBindTexture(GL_TEXTURE_1D, ctobj.texture);
        glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
        glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
        //glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
        glTexImage1D(GL_TEXTURE_1D, 0, GL_R32I, ctptr->facerightleftwithcolor.size(),
            0, GL_RED_INTEGER, GL_INT,ctptr->facerightleftwithcolor.data() );



        return;
    }

    void Scene::LoadLR(const MyLR& lr)
    {
        lrptr = &lr;
        SceneRenderLR temp;
        glGenVertexArrays(1, &temp.VAO);
        glGenBuffers(1, &temp.VBO);
        glGenBuffers(1, &temp.EBO);
        lrobj = temp;

        glActiveTexture(GL_TEXTURE1);
        glGenTextures(1, &lrobj.texture);
        glBindTexture(GL_TEXTURE_1D, lrobj.texture);
        glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
        glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
        //glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
        glTexImage1D(GL_TEXTURE_1D, 0, GL_R32I, lrptr->facetypeid.size(),
            0, GL_RED_INTEGER, GL_INT, lrptr->facetypeid.data());
        return;
    }

    void Scene::LoadLaceWire(const MyMesh& mesh, const MyCluster& clu)
    {
        const std::vector<LaceWires>& ref = clu.lacewires;

        glGenVertexArrays(1, &lineobj.VAO);
        glGenBuffers(1, &lineobj.VBO);
        glGenBuffers(1, &lineobj.EBO);
        uint i = 0;
        for (const auto& wireloop : ref) {
            for (const auto& wire : wireloop.wires) {
                for ( auto it = wire.cbegin(); it != wire.cend(); it++) {
                    uint idx = *it;
                    auto vh = mesh.vertex_handle(idx);
                    auto pnt = mesh.point(vh);
                    lineobj.geoinfo.emplace_back(pnt[0], pnt[1], pnt[2]);
                    lineobj.indices.push_back(i++);
                }
                //lineobj.geoinfo.pop_back();
                //lineobj.indices.pop_back();
                //i -= 1;
            }
            lineobj.indices.push_back(0xffffffff);
        }

    }

    void Scene::LoadLaceWire(const LaceWireGenerator& lwn,const MyMesh& mesh)
    {

        glGenVertexArrays(1, &lineobj.VAO);
        glGenBuffers(1, &lineobj.VBO);
        glGenBuffers(1, &lineobj.EBO);

        uint i = 0;
        for (const auto& wire : lwn.Ewires) {
            for (const auto& elem : wire.wire) {
                auto vh = mesh.vertex_handle(elem);
                auto pnt = mesh.point(vh);
                lineobj.geoinfo.emplace_back(pnt[0], pnt[1], pnt[2]);
                lineobj.indices.push_back(i++);
            }
            lineobj.indices.push_back(0xffffffff);
        }

    }

    void Scene::LoadLaceWire(const MyMesh& mesh, const NewCluster& nclu)
    {
        glGenVertexArrays(1, &lineobj.VAO);
        glGenBuffers(1, &lineobj.VBO);
        glGenBuffers(1, &lineobj.EBO);

        uint i = 0;
        for (const auto& wire : nclu.gloewires) {
            for (const auto& ver : wire) {
                auto vh = mesh.vertex_handle(ver);
                auto pnt = mesh.point(vh);
                lineobj.geoinfo.emplace_back(pnt[0], pnt[1], pnt[2]);
                lineobj.indices.push_back(i++);
            }
            lineobj.indices.push_back(0xffffffff);
        }


    }

    void Scene::LoadLacePoint(const MyMesh& mesh, const MyCluster& clu)
    {
        const std::vector<VertexSet> vsetref = clu.vertexsets;
        glGenVertexArrays(1, &pntobj.VAO);
        glGenBuffers(1, &pntobj.VBO);

        for (const auto& elems : vsetref) {
            //int idx = elem;
            //auto vh = mesh.vertex_handle(idx);
            //auto pnt = mesh.point(vh);
            //pntobj.geoinfo.emplace_back(pnt[0], pnt[1], pnt[2]);
            if (elems.boundset.size() != 0) {
                for (const auto& elem : elems.cornerset) {
                        int idx = elem;
                        auto vh = mesh.vertex_handle(idx);
                        auto pnt = mesh.point(vh);
                        pntobj.geoinfo.emplace_back(pnt[0], pnt[1], pnt[2]);
                }
                break;
            }
        }

    }

    void Scene::LoadInternalWire(const MyMesh& mesh, const LaceWireGenerator& lwn)
    {
        int idxcnt = 0;
        for (const auto& meshlet : lwn.meshlets) {
            for (const auto& elem : meshlet.interwire.wire) {
                if (elem == 0xFFFFFFFF) {
                    interobj.indices.push_back(0xffffffff);
                    continue;
                }
                else {
                    auto vh = mesh.vertex_handle(elem);
                    auto pnt = mesh.point(vh);
                    interobj.geoinfo.emplace_back(pnt[0], pnt[1], pnt[2]);
                    interobj.indices.push_back(idxcnt++);
                }
            }
            interobj.indices.push_back(0xffffffff);
        }
        glGenVertexArrays(1, &interobj.VAO);
        glGenBuffers(1, &interobj.VBO);
        glGenBuffers(1, &interobj.EBO);

    }

    void Scene::LoadInternalWire(const MyMesh& mesh, const NewLWGenerator& nlwn)
    {
        int idxcnt = 0;
        for (const auto& meshlet : nlwn.targets) {
            const auto& wire = meshlet.InternalLW;
            for (int i = 0; i < wire.vertex.size(); ++i) {
                if (wire.vertex[i] == EMPTYWIRE)
                    continue;
                auto vh = mesh.vertex_handle(wire.vertex[i]);
                auto pnt = mesh.point(vh);
                interobj.geoinfo.emplace_back(pnt[0], pnt[1], pnt[2]);
                interobj.indices.push_back(idxcnt++);
                if (wire.left[i] == WIREEND)
                    interobj.indices.push_back(0xFFFFFFFF);
            }
        }
        glGenVertexArrays(1, &interobj.VAO);
        glGenBuffers(1, &interobj.VBO);
        glGenBuffers(1, &interobj.EBO);

    }

    void Scene::LoadTSMeshlet(MyMesh mesh, Meshlets meshlets)
    {
        TS_MeshletLoad(tsobj.tsmeshlets, mesh, meshlets, tsobj.tsgeoinfo);

        //布局有问题
        tsobj.tsmeshlets[0].color = vec4(1.0f, 1.0f, 0.0f,1.0f);//test setting
        glGenBuffers(1, &tsobj.tsgeo);
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, tsobj.tsgeo);
        glBufferData(GL_SHADER_STORAGE_BUFFER, tsobj.tsgeoinfo.size() * sizeof(vec4), tsobj.tsgeoinfo.data(), GL_STATIC_DRAW);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, tsobj.tsgeo);
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);

        std::cout << int(&tsobj.tsgeoinfo[0]) - int(&tsobj.tsgeoinfo[1]) << std::endl;
        std::cout << sizeof(vec4) << std::endl;

        glGenBuffers(1, &tsobj.tsmeshdes);
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, tsobj.tsmeshdes);
        glBufferData(GL_SHADER_STORAGE_BUFFER, tsobj.tsmeshlets.size() * sizeof(TS_meshlet), tsobj.tsmeshlets.data(), GL_STATIC_DRAW);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 2, tsobj.tsmeshdes);
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);

        std::cout << int(&tsobj.tsmeshlets[0]) - int(&tsobj.tsmeshlets[1]) << std::endl;
        std::cout << sizeof(TS_meshlet) << std::endl;
    }

    void Scene::LoadIXMeshlet(MyMesh mesh, Meshlets meshlets)
    {
        IX_meshletLoad(ixobj.ixmeshlets, mesh, meshlets, ixobj.ixgeoinfo, ixobj.ixveridx, ixobj.ixprimidx);
        ixobj.ixmeshlets[1].color = vec4(1.0f, 1.0f, 0.0f, 1.0f);

        //几何信息
        glGenBuffers(1, &ixobj.IXGEO);
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, ixobj.IXGEO);
        glBufferData(GL_SHADER_STORAGE_BUFFER, ixobj.ixgeoinfo.size() * sizeof(vec4), ixobj.ixgeoinfo.data(), GL_STATIC_DRAW);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, ixobj.IXGEO);
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);
        //VER info
        glGenBuffers(1, &ixobj.IXVER);
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, ixobj.IXVER);
        glBufferData(GL_SHADER_STORAGE_BUFFER, ixobj.ixveridx.size() * sizeof(int), ixobj.ixveridx.data(), GL_STATIC_DRAW);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 2, ixobj.IXVER);
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);
        //PRIM info
        glGenBuffers(1, &ixobj.IXPRIM);
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, ixobj.IXPRIM);
        glBufferData(GL_SHADER_STORAGE_BUFFER, ixobj.ixprimidx.size() * sizeof(int), ixobj.ixprimidx.data(), GL_STATIC_DRAW);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 3, ixobj.IXPRIM);
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);
        //meshlet info
        glGenBuffers(1, &ixobj.IXMESHLET);
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, ixobj.IXMESHLET);
        glBufferData(GL_SHADER_STORAGE_BUFFER, ixobj.ixmeshlets.size() * sizeof(IX_meshlet), ixobj.ixmeshlets.data(), GL_STATIC_DRAW);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 4, ixobj.IXMESHLET);
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);



    }

    void Scene::LoadSCMeshlet(MyMesh mesh, Meshlets meshlets)
    {
        SC_meshletLoad(scobj.scmeshlets, mesh, meshlets, scobj.scgeoinfo, scobj.scprimidx);
        //几何信息
        glGenBuffers(1, &scobj.SCGEO);
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, scobj.SCGEO);
        glBufferData(GL_SHADER_STORAGE_BUFFER, scobj.scgeoinfo.size() * sizeof(vec4), scobj.scgeoinfo.data(), GL_STATIC_DRAW);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, scobj.SCGEO);
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);
        //PRIM info
        glGenBuffers(1, &scobj.SCPRIM);
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, scobj.SCPRIM);
        glBufferData(GL_SHADER_STORAGE_BUFFER, scobj.scprimidx.size() * sizeof(int), scobj.scprimidx.data(), GL_STATIC_DRAW);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 2, scobj.SCPRIM);
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);
        //meshlet info
        glGenBuffers(1, &scobj.SCMESHLET);
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, scobj.SCMESHLET);
        glBufferData(GL_SHADER_STORAGE_BUFFER, scobj.scmeshlets.size() * sizeof(SC_meshlet), scobj.scmeshlets.data(), GL_STATIC_DRAW);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 3, scobj.SCMESHLET);
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);

    }

    void Scene::LoadSimpleWireMeshlet(const LaceWireGenerator& lwn)
    {
        swobj.swgeoinfo = &lwn.geoinfo;
        swobj.swprimidx = &lwn.priminfo;
        swobj.simplemeshlets = &lwn.simplemeshlets;

        glGenBuffers(1, &swobj.SWGEO);
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, swobj.SWGEO);
        glBufferData(GL_SHADER_STORAGE_BUFFER, swobj.swgeoinfo->size() * sizeof(vec4), swobj.swgeoinfo->data(), GL_STATIC_DRAW);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, swobj.SWGEO);
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);
        //PRIM info
        glGenBuffers(1, &swobj.SWPRIM);
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, swobj.SWPRIM);
        glBufferData(GL_SHADER_STORAGE_BUFFER, swobj.swprimidx->size() * sizeof(unsigned char), swobj.swprimidx->data(), GL_STATIC_DRAW);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 2, swobj.SWPRIM);
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);
        //meshlet info
        glGenBuffers(1, &swobj.SWMESHLET);
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, swobj.SWMESHLET);
        glBufferData(GL_SHADER_STORAGE_BUFFER, swobj.simplemeshlets->size() * sizeof(Simple_meshlet), swobj.simplemeshlets->data(), GL_STATIC_DRAW);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 3, swobj.SWMESHLET);
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);

    }

    void Scene::LoadSimpleWireMeshlet(const NewLWGenerator& nlwn)
    {
        swobj.swgeoinfo = &nlwn.geoinfo;
        swobj.swprimidx = &nlwn.priminfo;
        swobj.simplemeshlets = &nlwn.simplemeshlets;

        glGenBuffers(1, &swobj.SWGEO);
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, swobj.SWGEO);
        glBufferData(GL_SHADER_STORAGE_BUFFER, swobj.swgeoinfo->size() * sizeof(vec4), swobj.swgeoinfo->data(), GL_STATIC_DRAW);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, swobj.SWGEO);
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);
        //PRIM info
        glGenBuffers(1, &swobj.SWPRIM);
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, swobj.SWPRIM);
        glBufferData(GL_SHADER_STORAGE_BUFFER, swobj.swprimidx->size() * sizeof(unsigned char), swobj.swprimidx->data(), GL_STATIC_DRAW);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 2, swobj.SWPRIM);
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);
        //meshlet info
        glGenBuffers(1, &swobj.SWMESHLET);
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, swobj.SWMESHLET);
        glBufferData(GL_SHADER_STORAGE_BUFFER, swobj.simplemeshlets->size() * sizeof(Simple_meshlet), swobj.simplemeshlets->data(), GL_STATIC_DRAW);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 3, swobj.SWMESHLET);
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);
    }


    void Scene::renderCT() {

        if (ctobj.needrender) {
            glBindVertexArray(ctobj.VAO);

            glBindBuffer(GL_ARRAY_BUFFER, ctobj.vertexVBO);
            glBufferData(GL_ARRAY_BUFFER, 3 * sizeof(float) * ctptr->vertices.size(), ctptr->vertices.data(), GL_DYNAMIC_DRAW);
            glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
            glEnableVertexAttribArray(0);

            glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ctobj.EBO);
            glBufferData(GL_ELEMENT_ARRAY_BUFFER,ctptr->Vtable.size()*sizeof(int), ctptr->Vtable.data(), GL_DYNAMIC_DRAW);

            //glBindTexture(GL_TEXTURE_1D, ctobj.texture);


            glDrawElements(GL_TRIANGLES, ctptr->Vtable.size(), GL_UNSIGNED_INT, 0);

        }
    }

    void Scene::renderML() {

    }

    void Scene::renderLR() {
        glBindVertexArray(lrobj.VAO);
        
        glBindBuffer(GL_ARRAY_BUFFER, lrobj.VBO);
        glBufferData(GL_ARRAY_BUFFER, 3 * sizeof(float) * lrptr->vertices.size(), lrptr->vertices.data(), GL_DYNAMIC_DRAW);
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
        glEnableVertexAttribArray(0);

        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, lrobj.EBO);
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, lrptr->EBO_check.size() * sizeof(int), lrptr->EBO_check.data(), GL_DYNAMIC_DRAW);

        glDrawElements(GL_TRIANGLES, lrptr->EBO_check.size(), GL_UNSIGNED_INT, 0);
    }

    void Scene::renderWireLine()
    {
        glBindVertexArray(lineobj.VAO);

        glBindBuffer(GL_ARRAY_BUFFER, lineobj.VBO);
        glBufferData(GL_ARRAY_BUFFER, 3 * sizeof(float) * lineobj.geoinfo.size(), lineobj.geoinfo.data(), GL_DYNAMIC_DRAW);
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
        glEnableVertexAttribArray(0);
        glEnable(GL_PRIMITIVE_RESTART_FIXED_INDEX);
        glLineWidth(5.0f);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, lineobj.EBO);
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, lineobj.indices.size() * sizeof(uint), lineobj.indices.data(), GL_DYNAMIC_DRAW);
        glDrawElements(GL_LINE_STRIP, lineobj.indices.size(), GL_UNSIGNED_INT, 0);
    }

    void Scene::renderInterWire()
    {
        glBindVertexArray(interobj.VAO);

        glBindBuffer(GL_ARRAY_BUFFER, interobj.VBO);
        glBufferData(GL_ARRAY_BUFFER, 3 * sizeof(float) * interobj.geoinfo.size(), interobj.geoinfo.data(), GL_DYNAMIC_DRAW);
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
        glEnableVertexAttribArray(0);
        glEnable(GL_PRIMITIVE_RESTART_FIXED_INDEX);
        glLineWidth(5.0f);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, interobj.EBO);
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, interobj.indices.size() * sizeof(uint), interobj.indices.data(), GL_DYNAMIC_DRAW);
        glDrawElements(GL_LINE_STRIP, interobj.indices.size(), GL_UNSIGNED_INT, 0);
    }

    void Scene::renderWirePnt()
    {
        glBindVertexArray(pntobj.VAO);
        glBindBuffer(GL_ARRAY_BUFFER, pntobj.VBO);
        glBufferData(GL_ARRAY_BUFFER, 3 * sizeof(float) * pntobj.geoinfo.size(), pntobj.geoinfo.data(), GL_DYNAMIC_DRAW);
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
        glEnableVertexAttribArray(0);
        glPointSize(10.0f);
        glDrawArrays(GL_POINTS, 0, pntobj.geoinfo.size());
    }

    void Scene::renderTSML()
    {
        glDrawMeshTasksNV(0, tsobj.tsmeshlets.size());
    }

    void Scene::renderIXML()
    {
        glDrawMeshTasksNV(0, ixobj.ixmeshlets.size());
    }

    void Scene::renderSCML()
    {
        glDrawMeshTasksNV(0, scobj.scmeshlets.size());
    }

    void Scene::renderSWML()
    {

        glDrawMeshTasksNV(0, swobj.simplemeshlets->size());

    }



    void Scene::renderCTring()
    {
        if (ctobj.ringneedrender) {
            glBindVertexArray(ctobj.VAO);

            glBindBuffer(GL_ARRAY_BUFFER, ctobj.vertexVBO);

            std::vector<int> temp = ctptr->ringvertex;
            temp.push_back(temp.front());
            glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ctobj.lineEBO);
            glBufferData(GL_ELEMENT_ARRAY_BUFFER, temp.size() * sizeof(int), temp.data(), GL_DYNAMIC_DRAW);
            glLineWidth(1.0f);
            glDrawElements(GL_LINE_STRIP, temp.size(), GL_UNSIGNED_INT, 0);

        }
    }

}
