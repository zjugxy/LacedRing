#include "scene.h"
#include "viewer.h"


namespace glfwviewer {
	void getGLError()
	{

	}


	Viewer::Viewer() :
		m_lightDirection(0, 0, -1, 0),
		m_rotateLightDir(0),
		m_scene(NULL)
	{

	}
	Viewer::~Viewer()
	{

	}
	void Viewer::set(Scene* pScene)
	{
		m_scene = pScene;
		vec3 center = scene()->bdSphere()->center();
		float radius = scene()->bdSphere()->radius();
		float viewdist = radius * 3;
		vec3 viewdir = vec3(0, 0, -1);
		camera()->setTarget(viewdir);
		camera()->setCenter(center - viewdir * viewdist);
		camera()->setUp(vec3(0, 1, 0));
		camera()->frustum()->setNearplane(radius / 10);
		camera()->frustum()->setFarplane(radius * 10);
		camera()->mover()->setStepsize(radius / 20);
		m_light.lightpos = camera()->center()+m_camera.up()*radius;
	}


	void Viewer::initGLFW() {

		GLFWwindow* window = NULL;
		int width, height;

		if (!glfwInit())
		{
			fprintf(stderr, "Failed to initialize GLFW\n");
			exit(EXIT_FAILURE);
		}

		glfwWindowHint(GLFW_DEPTH_BITS, 16);
		glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
		glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
		glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);


		window = glfwCreateWindow(1000, 1000, "ModelViewer", NULL, NULL);
		if (!window)
		{
			fprintf(stderr, "Failed to open GLFW window\n");
			glfwTerminate();
			exit(EXIT_FAILURE);
		}

		glfwMakeContextCurrent(window);
		// glad: load all OpenGL function pointers
		// ---------------------------------------
		if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
		{
			std::cout << "Failed to initialize GLAD" << std::endl;
			return;
		}
		glfwSwapInterval(1);
		m_window = window;

		static Viewer* v = this;
		//call back function of glfw
		glfwSetFramebufferSizeCallback(m_window, [](GLFWwindow* window, int width, int height) {
			v->reshape(window, width, height);
		});
		//get the mouse botton
		glfwSetMouseButtonCallback(m_window, [](GLFWwindow* window, int button, int action, int mods) {
			v->mouse_button_callback(window, button, action, mods);
		});
		//if pressed,then this func work
		glfwSetCursorPosCallback(m_window, [](GLFWwindow* window, double x, double y) {
			v->mousecallback(window, x, y);
		});
		glEnable(GL_DEPTH_TEST);
	}

	void Viewer::mouse_button_callback(GLFWwindow* window, int button, int action, int mods) {
		if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS) {
			sign_view = true;
			firstMouse = true;
		}
		if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_RELEASE) {
			sign_view = false;
			firstMouse = false;
		}
	}

	GLFWwindow* Viewer::MYwindow() const
	{
		return m_window;
	}

	void Viewer::mousecallback(GLFWwindow* window, double x, double y) {
		if (sign_view) {
			if (firstMouse) {
				Viewport* vp = m_camera.viewport();
				m_camera.tracker()->begin(x, y);
				firstMouse = false;
				return;
			}
			else
			{
				m_camera.tracker()->track(x, y);
				return;
			}
		}
	}

	void Viewer::reshape(GLFWwindow* window, int width, int height)
	{
		glViewport(0, 0, width, height);
		camera()->frustum()->setRatio(width / (float)height);
		camera()->viewport()->set(0, 0, width, height);
	}


	void Viewer::renderScene()
	{

		// draw the scene
		//model = glm::scale(model, vec3(uidata.myfloat, uidata.myfloat, uidata.myfloat));
		//model = glm::translate(model, vec3(0.2, 0.2, 0.0));
		m_shader.setMat4("projection", m_camera.projMatrix());
		m_shader.setMat4("view", m_camera.viewMatrix());
		//m_scene->draw()只渲染demo
		//if (uidata.drawtriangle)
		//	m_scene->draw();
		m_shader.setMat4("model", glm::mat4(1.0f));
		//m_scene->renderobjs();
		m_scene->renderCT();
	}

	void Viewer::processinput() {
		if (glfwGetKey(m_window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
			glfwSetWindowShouldClose(m_window, true);
		if (glfwGetKey(m_window, GLFW_KEY_W) == GLFW_PRESS)
			m_camera.mover()->moveForward();
		if (glfwGetKey(m_window, GLFW_KEY_S) == GLFW_PRESS)
			m_camera.mover()->moveBackward();
		if (glfwGetKey(m_window, GLFW_KEY_A) == GLFW_PRESS)
			m_camera.mover()->moveLeft();
		if (glfwGetKey(m_window, GLFW_KEY_D) == GLFW_PRESS)
			m_camera.mover()->moveRight();
	}

	void Viewer::display()
	{
		glClearColor(0.2f, 0.3f, 0.3f, 1.0f);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		m_shader.use();
		//m_shader.setBool("setcolor", uidata.setcolor);
		//m_shader.setVec3("color", uidata.mycolor[0], uidata.mycolor[1], uidata.mycolor[2]);

		m_shader.setVec3("viewPos", m_camera.center());
		m_shader.setVec3("lightPos", m_light.lightcolor);
		m_shader.setVec3("lightColor", m_light.lightcolor);

		renderScene();

		//修改
		//if(m_scene->)
	}



	void Viewer::printhelp(char* text[])
	{

	}

	void Viewer::Renderct(bool enabletyperender)
	{

		glClearColor(0.2f, 0.3f, 0.3f, 1.0f);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		m_shader.use();

		m_shader.setVec3("viewPos", m_camera.center());
		m_shader.setVec3("lightPos", m_light.lightpos);
		m_shader.setVec3("lightColor", m_light.lightcolor);

		glm::mat4 model = glm::mat4(1.0f);
		//model = glm::translate(model, vec3(0.2, 0.2, 0.0));
		m_shader.setMat4("projection", m_camera.projMatrix());
		m_shader.setMat4("view", m_camera.viewMatrix());
		m_shader.setMat4("model", model);
		//m_scene->draw()只渲染demo

		//unsigned int uniformfaceblock = glGetUniformBlockIndex(m_shader.ID, "faceidtocolor");
		//glUniformBlockBinding(m_shader.ID, uniformfaceblock, 0);
		//if(enabletyperender==false)
		//	m_shader.setArrayInt("faceidtocolor", m_scene->ctptr->facerightleftwithcolor);
		//else
		//	m_shader.setArrayInt("faceidtocolor", m_scene->ctptr->RenderTriType);

		m_shader.setInt("tex0", 0);
		m_shader.setBool("UseCT", true);

		m_scene->renderCT();

	}
	//运气很棒，可以使用一样的shader
	void Viewer::RenderuseLR() {
		glClearColor(0.2f, 0.3f, 0.3f, 1.0f);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		m_shader.use();

		m_shader.setVec3("viewPos", m_camera.center());
		m_shader.setVec3("lightPos", m_light.lightpos);
		m_shader.setVec3("lightColor", m_light.lightcolor);

		glm::mat4 model = glm::mat4(1.0f);
		//model = glm::translate(model, vec3(0.2, 0.2, 0.0));
		m_shader.setMat4("projection", m_camera.projMatrix());
		m_shader.setMat4("view", m_camera.viewMatrix());
		m_shader.setMat4("model", model);
		//m_shader.setArrayInt("faceidtocolor", m_scene->lrptr->facetypeid);
		m_shader.setBool("UseCT", false);
		m_shader.setInt("tex1", 1);

		m_scene->renderLR();

	}

	void Viewer::RenderML()
	{
		glClearColor(0.2f, 0.3f, 0.3f, 1.0f);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		m_shader.use();

		m_shader.setVec3("viewPos", m_camera.center());
		m_shader.setVec3("lightPos", m_light.lightpos);
		m_shader.setVec3("lightColor", m_light.lightcolor);

		glm::mat4 model = glm::mat4(1.0f);
		//model = glm::translate(model, vec3(0.2, 0.2, 0.0));
		m_shader.setMat4("projection", m_camera.projMatrix());
		m_shader.setMat4("view", m_camera.viewMatrix());
		m_shader.setMat4("model", model);
		//m_shader.setArrayInt("faceidtocolor", m_scene->lrptr->facetypeid);
		m_shader.setBool("UseCT", false);
		m_shader.setInt("tex2", 2);

		m_shader.setInt("n_face", m_scene->mlobj.meshptr->n_faces());

		m_scene->renderML();
	}

	void Viewer::Renderring()
	{
		ring_shader.use();
		glm::mat4 model = glm::mat4(1.0f);
		//model = glm::translate(model, vec3(0.2, 0.2, 0.0));
		ring_shader.setMat4("projection", m_camera.projMatrix());
		ring_shader.setMat4("view", m_camera.viewMatrix());
		ring_shader.setMat4("model", model);
		m_scene->renderCTring();
	}

	

} //namespace glfwviewer