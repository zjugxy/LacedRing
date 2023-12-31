#include<memory>
#include "scene.h"

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
        //����mmbox mmsphere
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

    //void Scene::renderobjs() {
    //    for (auto& obj : renlist)
    //        if (obj.needrender) {
    //            obj.ptr->Generate_Rendata();
    //            glBindVertexArray(obj.VAO);

    //            glBindBuffer(GL_ARRAY_BUFFER, obj.vertexVBO);
    //            glBufferData(GL_ARRAY_BUFFER, 3 * sizeof(float) * obj.ptr->renderdata.points.size(), obj.ptr->renderdata.points.data(), GL_DYNAMIC_DRAW);
    //            glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
    //            glEnableVertexAttribArray(0);


    //            glBindBuffer(GL_ARRAY_BUFFER, obj.normalVBO);
    //            glBufferData(GL_ARRAY_BUFFER, 3 * sizeof(float) * obj.ptr->renderdata.normals.size(), obj.ptr->renderdata.normals.data(), GL_DYNAMIC_DRAW);
    //            glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
    //            glEnableVertexAttribArray(1);

    //            glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, obj.EBO);
    //            glBufferData(GL_ELEMENT_ARRAY_BUFFER, obj.ptr->renderdata.faces.size() * sizeof(int), obj.ptr->renderdata.faces.data(), GL_DYNAMIC_DRAW);
    //            //����
    //            glDrawElements(GL_TRIANGLES, obj.ptr->renderdata.faces.size(), GL_UNSIGNED_INT, 0);
    //        }

    //}

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
