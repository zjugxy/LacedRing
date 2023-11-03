#include<memory>
#include "scene.h"
#include<random>
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

    void Scene::LoadMeshlet(MyMesh mesh, std::vector<std::vector<int>> meshlets)
    {



        for (auto& vertex : mesh.vertices()) {
            auto pnt = mesh.point(vertex);
            mlobj.vertices.emplace_back(pnt[0], pnt[1], pnt[2]);
        }
        for (auto& face : mesh.faces()) {
            for (const auto& ver : face.vertices()) {
                mlobj.indexs.push_back(ver.idx());
            }
        }

        mlobj.meshptr = &mesh;
        glGenVertexArrays(1, &mlobj.VAO);
        glGenBuffers(1, &mlobj.VBO);
        glGenBuffers(1, &mlobj.EBO);

        glActiveTexture(GL_TEXTURE2);
        glGenTextures(1, &mlobj.texture);
        glBindTexture(GL_TEXTURE_1D, mlobj.texture);
        glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
        glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);



        // 绑定颜色数组到纹理对象
        glTexImage1D(GL_TEXTURE_1D, 0, GL_RGB, color.size(), 0, GL_RGB, GL_UNSIGNED_BYTE, color.data());

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
        glBindVertexArray(mlobj.VAO);

        glBindBuffer(GL_ARRAY_BUFFER, mlobj.VBO);
        glBufferData(GL_ARRAY_BUFFER, 3 * sizeof(float) * mlobj.vertices.size(), mlobj.vertices.data(), GL_DYNAMIC_DRAW);
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
        glEnableVertexAttribArray(0);

        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, mlobj.EBO);
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, mlobj.indexs.size()*sizeof(uint), mlobj.indexs.data(), GL_DYNAMIC_DRAW);

        glDrawElements(GL_TRIANGLES, mlobj.indexs.size(), GL_UNSIGNED_INT, 0);
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
