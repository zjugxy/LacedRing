#ifndef _______viewer_camera_h___________
#define _______viewer_camera_h___________

#include"viewermath.h"

namespace glfwviewer {

	class Viewport;
	class Frustum;
	class Camera;

	/** */
	class Viewport {
	public:
		Viewport();
		void set(int x, int y, int width, int height);
		int x() const { return m_x; }
		int y() const { return m_y; }
		int width() const { return m_width; }
		int height() const { return m_height; }
		int halfwidth() const { return m_width / 2; }
		int halfheight() const { return m_height / 2; }

	protected:
		int m_x;
		int m_y;
		int m_width;
		int m_height;
	};

	class Frustum {
	public:
		Frustum();
		/** */
		float fovy() const { return m_fovy; }
		float ratio() const { return m_ratio; }
		float nearplane() const { return m_near; }
		float farplane() const { return m_far; }
		void  setFovy(float val) { m_fovy = val; }
		void  setRatio(float val) { m_ratio = val; }
		void  setNearplane(float val) { m_near = val; }
		void  setFarplane(float val) { m_far = val; }
		void  set(float vfovy, float vratio, float vnear, float vfar);
	protected:
		float m_fovy;
		float m_ratio;
		float m_near;
		float m_far;
		vec3  m_points[8];

	};

	/** */
	class CameraTracker {
	public:
		CameraTracker();
		void setCamera(Camera* cam) { m_camera = cam; }
		void begin(int x, int y);
		void track(int x, int y);
		void track(int x, int y, vec3& v);
	protected:
		Camera* m_camera;
		vec2 m_prev_point;
		float m_sensitivity;
	};

	/** */
	class CameraMover {
	public:
		CameraMover();
		void    setCamera(Camera* cam) { m_camera = cam; }
		Camera* camera() { return m_camera; }
		void    setStepsize(float stepsize) { m_stepsize = stepsize; }
		float   stepsize() const { return m_stepsize; }
		void    moveForward();
		void    moveBackward();
		void    moveLeft();
		void    moveRight();
	protected:
		Camera* m_camera;
		float m_stepsize;
	};

	/** */
	class Camera {
	public:
		Camera();
		~Camera();
		/** Rotates camera vectors by angle, given axis */
		void rotate(float angle, const vec3& axis);
		/** Returns the view matrix of the camera */
		mat4 viewMatrix();
		/** Returns the projection matrix defined by the current Frustum */
		mat4 projMatrix();

		vec3 target() const { return m_target; }
		vec3 center() const { return m_center; }
		vec3 up() const { return m_up; }
		vec3 right() const;
		vec3 left() const;

		void setTarget(const vec3& target) { m_target = target; }
		void setCenter(const vec3& center) { m_center = center; }
		void setUp(const vec3& up) { m_up = up; }

		Frustum* frustum() { return &m_frustum; }
		Viewport* viewport() { return &m_viewport; }
		CameraTracker* tracker() { return &m_tracker; }
		CameraMover* mover() { return &m_mover; }

	private:
		vec3 m_center;   // camera center position
		vec3 m_target;   // viewing direction
		vec3 m_up;       // up direction
		Viewport m_viewport;
		Frustum  m_frustum;
		CameraTracker m_tracker;
		CameraMover   m_mover;
	};

}

#endif
