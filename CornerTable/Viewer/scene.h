#pragma once
#ifndef ______scene_h_____________________
#define ______scene_h_____________________
#include <glad/glad.h>
#include<GLFW/glfw3.h>
#include "viewermath.h"
#include"../CornerTable.h"
#include"../MyLR.h"
#include"../Meshlet.h"
namespace glfwviewer {

	class BoundingSphere
	{
		vec3 m_center;
		float m_radius;
	public:
		BoundingSphere(float x = 0, float y = 0, float z = 0, float radius = 1)
			: m_center(x, y, z), m_radius(radius)
		{
		}
		void center(vec3 cen) { m_center = cen; }
		void radius(float r) { m_radius = r; }
		vec3 center() const { return m_center; }
		float radius() const { return m_radius; }
	};

	class BoundingBox
	{
		vec3 m_min;
		vec3 m_max;
	public:
		BoundingBox()
			: m_min{ -1,-1,-1 }, m_max(1, 1, 1)
		{
		}
		BoundingBox(vec3 min, vec3 max)
			: m_min(min), m_max(max)
		{
		}
		vec3 mymin() const { return m_min; }
		vec3 mymax() const { return m_max; }
		void mymin(vec3 p) { m_min = p; }
		void mymax(vec3 p) { m_max = p; }
	};

	struct SceneRenderobj
	{
		GLuint normalVBO, vertexVBO;
		GLuint EBO;
		GLuint VAO;
		bool needrender = true;
	};

	struct SceneRenderCT {
		GLuint vertexVBO;
		GLuint EBO;
		GLuint VAO;
		GLuint facecolorVBO;//uniform?
		bool needrender = true;
		unsigned int texture;
		GLuint lineEBO;
		bool ringneedrender = true;

	};

	struct SceneRenderLR {
		GLuint VAO;
		GLuint VBO;
		GLuint EBO;
		unsigned int texture;


	};

	struct SceneRenderMeshlet {
		MyMesh* meshptr;
		std::vector<vec3> vertices;
		std::vector<uint> indexs;
		GLuint VAO;
		GLuint VBO;
		GLuint EBO;
		unsigned int texture;
	};

	struct SceneRenderTS {
		std::vector<TS_meshlet> tsmeshlets;
		std::vector<vec4> tsgeoinfo;
		GLuint tsgeo;
		GLuint tsmeshdes;
	};


	struct SceneRenderIX {
		std::vector<IX_meshlet> ixmeshlets;
		std::vector<vec4> ixgeoinfo;
		std::vector<int> ixveridx;
		std::vector<int> ixprimidx;
		GLuint IXGEO;
		GLuint IXVER;
		GLuint IXPRIM;
		GLuint IXMESHLET;
	};

	struct SceneRenderSC {
		std::vector<SC_meshlet> scmeshlets;
		std::vector<vec4> scgeoinfo;
		std::vector<int> scprimidx;
		GLuint SCGEO;
		GLuint SCPRIM;
		GLuint SCMESHLET;
	};


	class Scene
	{
	public:
		BoundingBox m_bbox;
		BoundingSphere m_bsphere;
		GLuint VAO;
		SceneRenderCT ctobj;
		const CornerTable* ctptr;
		const MyLR* lrptr;
		SceneRenderLR lrobj;
		SceneRenderMeshlet mlobj;
		std::vector<SceneRenderobj> renlist;

		SceneRenderTS tsobj;
		SceneRenderIX ixobj;
		SceneRenderSC scobj;

	public:
		Scene();
		~Scene();
		virtual void initialize();
		BoundingSphere* bdSphere() { return &m_bsphere; }
		BoundingBox* bdBox() { return &m_bbox; }

		void LoadCornerTable(const CornerTable& ct);
		void LoadLR(const MyLR& lr);

		void renderML();
		void renderCT();
		void renderCTring();
		void renderLR();

		//void renderobjs();


		//meshlet load
		void LoadTSMeshlet(MyMesh mesh, Meshlets meshlets);
		void LoadIXMeshlet(MyMesh mesh, Meshlets meshlets);
		void LoadSCMeshlet(MyMesh mesh, Meshlets meshlets);


		void renderTSML();
		void renderIXML();
		void renderSCML();
	};

}

#endif