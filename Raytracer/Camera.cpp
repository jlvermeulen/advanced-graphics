#define _USE_MATH_DEFINES

#include <Camera.h>
#include <math.h>

#define M_PI_180 M_PI/180.0f


//--------------------------------------------------------------------------------
Camera::Camera() :
  eye_(Vector3F(0, 0, 3)),
  focus_(Vector3F::Forward),
  right_(Vector3F::Right),
  up_(Vector3F::Up),
  Width(1280),
  Height(720),
  fovY_(40.0f),
  zNear_(0.1f),
  zFar_(20.0f),
  Aperture(3.0f),
  FocalDepth(4.0f)
{

}

//--------------------------------------------------------------------------------
Camera::Camera(Vector3F eye, Vector3F focus, Vector3F up, int width, int height, float fovY, float zNear, float zFar, float aperture, float focalDepth) :
  eye_(eye),
  focus_(focus),
  up_(up),
  Width(width),
  Height(height),
  fovY_(fovY),
  zNear_(zNear),
  zFar_(zFar),
  Aperture(aperture),
  FocalDepth(focalDepth)
{
  right_ = Vector3F::Cross(focus_, up_);
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
void Camera::Rotate(const Vector3F& axis, float angle)
{
  Quaternion rotation(
    axis.X * sin(angle / 2),
    axis.Y * sin(angle / 2),
    axis.Z * sin(angle / 2),
    cos(angle / 2)
  );

  Quaternion view(
    focus_.X,
    focus_.Y,
    focus_.Z,
    0
  );

  Quaternion result = (rotation * view) * Quaternion::conjugate(rotation);

  focus_.X = result.X;
  focus_.Y = result.Y;
  focus_.Z = result.Z;
}

//--------------------------------------------------------------------------------
void Camera::RotateX(float angle)
{
  //Rotate(Vector3F::Cross(focus_, up_), 0.01f * angle);

  focus_ = Vector3F::Normalise(focus_ * cosf(angle * M_PI_180) + up_ * sinf(angle * M_PI_180));
  up_ = Vector3F::Cross(focus_, right_) * -1;
}

//--------------------------------------------------------------------------------
void Camera::RotateY(float angle)
{
  //Rotate(Vector3F::Up, 0.01f * angle);

  focus_ = Vector3F::Normalise(focus_ * cosf(angle * M_PI_180) - right_ * sinf(angle * M_PI_180));
  right_ = Vector3F::Cross(focus_, up_);
}

//--------------------------------------------------------------------------------
void Camera::RotateZ(float angle)
{
  right_ = Vector3F::Normalise(right_ * cosf(angle * M_PI_180) + up_ * sinf(angle * M_PI_180));
  up_ = Vector3F::Cross(focus_, right_) * -1;
}