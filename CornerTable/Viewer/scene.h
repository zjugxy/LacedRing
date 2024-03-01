#pragma once
#ifndef ______scene_h_____________________
#define ______scene_h_____________________
#include <glad/glad.h>
#include<GLFW/glfw3.h>
#include "viewermath.h"
#include"../CornerTable.h"
#include"../MyLR.h"
#include"../Meshlet.h"
#include"../MyCluster.h"
#include"../NewCluster.h"
#include"../LaceWireGenerator.h"
#include"../NewLWGenerator.h"
#include"../MeshFile.h"


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

	struct SceneRenderWireline
	{
		GLuint VAO;
		GLuint VBO;
		GLuint EBO;
		std::vector<uint> indices;
		std::vector<vec3> geoinfo;
	};

	struct SceneRenderWirePoint {
		GLuint VAO;
		GLuint VBO;
		std::vector<vec3> geoinfo;
	};

	struct SceneRenderInterWire {
		GLuint VAO;
		GLuint VBO;
		GLuint EBO;
		std::vector<uint> indices;
		std::vector<vec3> geoinfo;
	};

	struct SceneRenderSimpleWire
	{
		const std::vector<Simple_meshlet>* simplemeshlets;
		const std::vector<vec4>* swgeoinfo;
		const std::vector<unsigned char>* swprimidx;
		GLuint SWGEO;
		GLuint SWPRIM;
		GLuint SWMESHLET;

	};

	struct SceneGPULW {
		const std::vector<uint>* DesLoc;
		const std::vector<uint>* Desinfo;

		const std::vector<uint>* newintercon;//inter: left,right,irregular
		const std::vector<uint>* newextercon;//exter: left,right
		const std::vector<float>* intergeo;
		const std::vector<float>* extergeo;

		GLuint LWdesloc, LWdesinfo, LWintercon, LWextercon, LWintergeo, LWextergeo;
	};

	struct SceneFinalGPULW {
		const std::vector<uint>* DesLoc;
		const std::vector<uint>* Desinfo;

		const std::vector<uint>* newintercon;//inter: left,right,irregular
		const std::vector<uint>* newextercon;//exter: left,right
		const std::vector<uint>* finalintergeo;
		const std::vector<uint>* finalextergeo;

		GLuint FinalLWdesloc, FinalLWdesinfo, FinalLWintercon, 
			FinalLWextercon, FinalLWintergeo, FinalLWextergeo;
	};


	struct SceneRenderNormal {
		std::vector<vec3> linepoints;
		GLuint VAO;
		GLuint VBO;
	};

	struct SceneRenderPoint {
		std::vector<vec3> points;
		GLuint VAO;
		GLuint VBO;
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

		SceneRenderWireline lineobj;
		SceneRenderWirePoint pntobj;
		SceneRenderInterWire interobj;
		SceneRenderSimpleWire swobj;
		SceneGPULW gpulwobj;
		SceneFinalGPULW finalgpulwobj;
		SceneRenderNormal normalobj;
		SceneRenderPoint pointobj;

	public:
		Scene();
		~Scene();
		virtual void initialize();
		BoundingSphere* bdSphere() { return &m_bsphere; }
		BoundingBox* bdBox() { return &m_bbox; }

		void LoadCornerTable(const CornerTable& ct);
		void LoadLR(const MyLR& lr);
		void LoadLaceWire(const MyMesh& mesh,const MyCluster& clu);
		void LoadLaceWire(const LaceWireGenerator& lwn,const MyMesh& mesh);
		void LoadLaceWire(const MyMesh& mesh, const NewCluster& nclu);

		void LoadLacePoint(const MyMesh& mesh, const MyCluster& clu);
		void LoadInternalWire(const MyMesh& mesh ,const LaceWireGenerator& lwn);
		void LoadInternalWire(const MyMesh& mesh, const NewLWGenerator& nlwn);
		void LoadNormalLine(const NewCluster& nclu);
		void LoadPoints(const MyMesh& mesh);
		void LoadLines(const MyMesh& mesh,Meshlets meshlet);


		void renderML();
		void renderCT();
		void renderCTring();
		void renderLR();

		void renderWireLine();
		void renderInterWire();
		void renderWirePnt();
		void renderNormalLine();
		void renderPoints();
		//void renderobjs();


		//meshlet load
		void LoadTSMeshlet(MyMesh mesh, Meshlets meshlets);
		void LoadIXMeshlet(MyMesh mesh, Meshlets meshlets);
		void LoadSCMeshlet(MyMesh mesh, Meshlets meshlets);
		void LoadSimpleWireMeshlet(const LaceWireGenerator& lwn);
		void LoadSimpleWireMeshlet(const NewLWGenerator& nlwn);
		void LoadGPULW(NewLWGenerator& nlwn);
		void LoadGPULW(MeshFile& meshfile);
		void LoadFinalLaceWire(NewLWGenerator& nlwn);

		void renderTSML();
		void renderIXML();
		void renderSCML();
		void renderSWML();
		void renderGPULW();
		void renderFinalGPULW();
	};

}

#endif