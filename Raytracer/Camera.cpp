#define _USE_MATH_DEFINES

#include <Camera.h>
#include <math.h>

Camera::Camera(Vector3D eye, Vector3D focus, Vector3D up, double aspectRatio) :
  eye_(eye),
  focus_(focus),
  up_(up)
{
  // Instantiate the projection matrix
  projection_ = Matrix4x4D::createPerspectiveFieldOfView(M_1_PI, aspectRatio, 1.0, 300.0);

  // Instantiate the view matrix
  updateViewMatrix();
}

void Camera::setEye(Vector3D eye)
{
  eye_ = eye;
  updateViewMatrix();
}

void Camera::setFocus(Vector3D focus)
{
  focus_ = focus;
  updateViewMatrix();
}

void Camera::setUp(Vector3D up)
{
  up_ = up;
  updateViewMatrix();
}

void Camera::updateViewMatrix()
{
  view_ = Matrix4x4D::createLookAt(eye_, focus_, up_);
}