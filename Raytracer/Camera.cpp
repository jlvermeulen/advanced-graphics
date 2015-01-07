#define _USE_MATH_DEFINES

#include <Camera.h>
#include <math.h>

#define M_PI_180 M_PI/180.0

//--------------------------------------------------------------------------------
Camera::Camera() :
  eye_(Vector3D(0, 0, 3)),
  focus_(Vector3D::Forward),
  right_(Vector3D::Right),
  up_(Vector3D::Up),
  Width(1280),
  Height(720),
  fovY_(40.0f),
  zNear_(0.1f),
  zFar_(20.0f)
{

}

//--------------------------------------------------------------------------------
Camera::Camera(Vector3D eye, Vector3D focus, Vector3D up, int width, int height, float fovY, float zNear, float zFar) :
  eye_(eye),
  focus_(focus),
  up_(up),
  Width(width),
  Height(height),
  fovY_(fovY),
  zNear_(zNear),
  zFar_(zFar)
{
  right_ = Vector3D::Cross(focus_, up_);
}

//--------------------------------------------------------------------------------
Camera::~Camera()
{

}

//--------------------------------------------------------------------------------
void Camera::MoveForward(const float& distance)
{
  eye_ += focus_ * distance;
}

//--------------------------------------------------------------------------------
void Camera::MoveRight(const float& distance)
{
  eye_ += right_ * distance;
}

//--------------------------------------------------------------------------------
void Camera::MoveUpward(const float& distance)
{
  eye_ += up_ * distance;
}

//--------------------------------------------------------------------------------
void Camera::RotateX(float angle)
{
  rotX_ += angle;

  focus_ = Vector3D::Normalise(focus_ * cos(angle * M_PI_180) + up_ * sin(angle * M_PI_180));
  up_ = Vector3D::Cross(focus_, right_) * -1;
}

//--------------------------------------------------------------------------------
void Camera::RotateY(float angle)
{
  rotY_ += angle;

  focus_ = Vector3D::Normalise(focus_ * cos(angle * M_PI_180) - right_ * sin(angle * M_PI_180));
  right_ = Vector3D::Cross(focus_, up_);
}

//--------------------------------------------------------------------------------
void Camera::RotateZ(float angle)
{
  rotZ_ += angle;

  right_ = Vector3D::Normalise(right_ * cos(angle * M_PI_180) + up_ * sin(angle * M_PI_180));
  up_ = Vector3D::Cross(focus_, right_) * -1;
}