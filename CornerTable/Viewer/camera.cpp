#include "camera.h"

/** */
namespace glfwviewer {

	/** */
	Viewport::Viewport() :
		m_x(0),
		m_y(0),
		m_width(0),
		m_height(0)
	{
	}

	/** */
	void Viewport::set(int x, int y, int width, int height)
	{
		m_x = x;
		m_y = y;
		m_width = width;
		m_height = height;
	}

	/** */
	Frustum::Frustum() :
		m_fovy(45.0),
		m_ratio(0.5),
		m_near(1.0),
		m_far(1000.0),
		m_points{}
	{

	}

	/** */
	void Frustum::set(float vfovy, float vratio, float vnear, float vfar)
	{
		m_fovy = vfovy;
		m_ratio = vratio;
		m_near = vnear;
		m_far = vfar;
	}


	/** */
	CameraTracker::CameraTracker() :
		m_camera(NULL),
		m_prev_point(0.0f),
		m_sensitivity(0.2f)
	{
	}

	/** */
	void CameraTracker::begin(int x, int y)
	{
		m_prev_point = vec2(x, y);
	}

	/** */
	void CameraTracker::track(int x, int y)
	{
		if (!m_camera) return;
		vec2 pnow = vec2(x, y);
		vec2 vrot = (pnow - m_prev_point) * m_sensitivity;
		m_camera->rotate(vrot.x, m_camera->up());
		m_camera->rotate(vrot.y, m_camera->right());
		m_prev_point = pnow;
	}

	void CameraTracker::track(int x, int y, vec3& v)
	{
		if (!m_camera) return;
		vec2 n_point = vec2(x, y);
		vec2 n_rotation = (n_point - m_prev_point) * m_sensitivity;
		v = glm::rotate(v, glm::radians(n_rotation.x), m_camera->up());
		v = glm::rotate(v, glm::radians(n_rotation.y), m_camera->right());
		v = glm::normalize(v);
		m_prev_point = n_point;
	}

	CameraMover::CameraMover() :
		m_camera(NULL),
		m_stepsize(1.0f)
	{
	}

	void CameraMover::moveForward()
	{
		if (!m_camera) return;
		camera()->setCenter(camera()->center()
			+ camera()->target() * stepsize());
	}

	void CameraMover::moveBackward()
	{
		if (!m_camera) return;
		camera()->setCenter(camera()->center()
			- camera()->target() * stepsize());
	}

	void CameraMover::moveLeft()
	{
		if (!m_camera) return;
		camera()->setCenter(camera()->center()
			+ camera()->left() * stepsize());
	}

	void CameraMover::moveRight()
	{
		if (!m_camera) return;
		camera()->setCenter(camera()->center()
			+ camera()->right() * stepsize());
	}

	/** */
	Camera::Camera() :
		m_center(0.0f),
		m_target(0.0f, 0.0f, 1.0f),
		m_up(0.0f, 1.0f, 0.0f)
	{
		m_tracker.setCamera(this);
		m_mover.setCamera(this);
	}

	Camera::~Camera()
	{
	}

	/** */
	void Camera::rotate(float angle, const vec3& axis)
	{
		m_target = glm::normalize(glm::rotate(
			m_target, glm::radians(angle), axis));
		m_up = glm::normalize(glm::rotate(
			m_up, glm::radians(angle), axis));
	}

	/** */
	mat4 Camera::viewMatrix()
	{
		mat4 vm = glm::lookAt(m_center, m_center + m_target, m_up);
		//printf("cam: %f %f %f cam.ratio:%f\n", m_center.x, m_center.y, m_center.z, this->frustum()->ratio());
		return vm;
	}

	/** */
	mat4 Camera::projMatrix()
	{
		return glm::perspective(m_frustum.fovy(), m_frustum.ratio(), m_frustum.nearplane(), m_frustum.farplane());
	}

	vec3 Camera::right() const
	{
		return glm::normalize(glm::cross(m_target, m_up));
	}

	vec3 Camera::left() const
	{
		return glm::normalize(glm::cross(m_up, m_target));
	}

}
