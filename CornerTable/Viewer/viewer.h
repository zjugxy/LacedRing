#pragma once
#include<glad/glad.h>
#include<GLFW/glfw3.h>
#include "viewermath.h"
#include "camera.h"
#include "../Shader.h"
#include"../MShader.h"
namespace glfwviewer {


	void getGLError();

	class Scene;

	struct UIdata
	{
		float myfloat = 1.0;
		bool drawtriangle = true;
		bool setcolor = false;
		float* mycolor = (float*)malloc(sizeof(float) * 3);
		bool showwindow = false;
	};

	struct lightdata {
		vec3 lightpos;//初始化为相机初始位置
		vec3 lightcolor = glm::vec3(1.0, 1.0, 1.0);
	};

	class Viewer
	{
	public:
		Viewer();
		~Viewer();
	public:
		virtual void set(Scene* scene);

		//virtual void guiInit();
		//virtual void guiQuit();
		//virtual void guiDisplay();
		virtual void renderScene();
		virtual void display();
		virtual void printhelp(char* text[]);
		virtual void Renderct(bool enabletyperender = false);
		virtual void RenderuseLR();
		virtual void RenderTSML();
		virtual void Renderring();
		virtual void RenderIXML();
		virtual void RenderSCML();
		virtual void RenderSWML();
		virtual void RenderFINALLWML();

		virtual void RenderWireLine();
		virtual void RenderWirePnt();
		virtual void RenderInterWire();
		//input func below
		virtual void initGLFW();
		virtual void processinput();
		virtual void mousecallback(GLFWwindow* window, double xposIn, double yposIn);
		virtual void reshape(GLFWwindow* window, int width, int height);
		virtual void mouse_button_callback(GLFWwindow* window, int button, int action, int mods);


		GLFWwindow* MYwindow() const;

		//add function
		virtual void setWindow(GLFWwindow* window) {
			m_window = window;
			int width, height;
			glfwGetFramebufferSize(window, &width, &height);
			m_camera.frustum()->setRatio(float(width) / height);
			m_camera.viewport()->set(0, 0, width, height);

		}
		virtual void setshader(const char* vertexPath, const char* fragmentPath, const char* geometryPath = nullptr) {
			m_shader = Shader(vertexPath, fragmentPath, geometryPath);
		}
		virtual void setlineshader(const char* vertexPath, const char* fragmentPath, const char* geometryPath = nullptr) {
			m_line_shader = Shader(vertexPath, fragmentPath, geometryPath);
		}
		virtual void setMeshshader(const char* mesh, const char* frag);
		virtual void setringshader(const char* vertexPath, const char* fragmentPath, const char* geometryPath = nullptr) {
			ring_shader = Shader(vertexPath, fragmentPath, geometryPath);
		}

	public:
		Camera* camera() { return &m_camera; }
		Scene* scene() { return m_scene; }
	protected:
		Scene* m_scene;
		Camera m_camera;
		vec4 m_lightDirection;
		int    m_rotateLightDir;

		GLFWwindow* m_window;
		Shader m_shader;
		Shader ring_shader;
		Shader LR_shader;
		MShader m_meshshader;
		Shader m_line_shader;
		bool firstMouse = true;
		bool sign_view = false;

	public:
		//UIdata uidata;
		lightdata m_light;
	};


}
