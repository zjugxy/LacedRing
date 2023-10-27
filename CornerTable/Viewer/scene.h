#pragma once
#ifndef ______scene_h_____________________
#define ______scene_h_____________________
#include <glad/glad.h>
#include<GLFW/glfw3.h>
#include "viewermath.h"
#include"../CornerTable.h"
#include"../MyLR.h"

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

		GLuint lineEBO;
		bool ringneedrender = true;

	};

	struct SceneRenderLR {
		GLuint VAO;
		GLuint VBO;
		GLuint EBO;
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
		
		std::vector<SceneRenderobj> renlist;

	public:
		Scene();
		~Scene();
		virtual void initialize();
		BoundingSphere* bdSphere() { return &m_bsphere; }
		BoundingBox* bdBox() { return &m_bbox; }

		void LoadCornerTable(const CornerTable& ct);
		void LoadLR(const MyLR& lr);

		void renderCT();
		void renderCTring();
		void renderLR();
		//void renderobjs();
	};

}

#endif