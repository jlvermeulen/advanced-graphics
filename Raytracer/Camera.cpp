#include <Camera.h>

#define PI 3.1415265359
#define PId 3.1415265359/180.0

//--------------------------------------------------------------------------------
Camera::Camera() :
  changed_(false),
  focus_(Vector3D::Forward)
{

}

//--------------------------------------------------------------------------------
Camera::~Camera()
{

}

//--------------------------------------------------------------------------------
void Camera::MoveForward(const float& distance)
{
  Vector3D viewDir = getViewDirection();

  Vector3D move;
  move.X = viewDir.X * -distance;
  move.Y = viewDir.Y * -distance;
  move.Z = viewDir.Z * -distance;

  eye_ += move;
}

//--------------------------------------------------------------------------------
void Camera::MoveSideways(const float& distance)
{
  Vector3D viewDir = getViewDirection();

  Vector3D move;
  move.X = viewDir.Z * -distance;
  move.Y = 0.0;
  move.Z = -viewDir.X * -distance;

  eye_ += move;
}

//--------------------------------------------------------------------------------
void Camera::Move(const Vector3D& direction)
{
  eye_ += direction;
}

//--------------------------------------------------------------------------------
void Camera::RotateX(float angle)
{
  changed_ = true;
  rotX_ += angle;
}

//--------------------------------------------------------------------------------
void Camera::RotateY(float angle)
{
  changed_ = true;
  rotY_ += angle;
}

//--------------------------------------------------------------------------------
void Camera::RotateZ(float angle)
{
  changed_ = true;
  rotZ_ += angle;
}

//--------------------------------------------------------------------------------
const Vector3D& Camera::getViewDirection()
{
  if (changed_)
  {
    Vector3D temp;

    // Rotate around Y-axis
    temp.X = cos((rotY_ + 90.0) * PId);
    temp.Z = -sin((rotY_ + 90.0) * PId);

    // Rotate around X-axis
    double cosX = cos(rotX_ * PId);
    view_.X = temp.X * cosX;
    view_.Y = sin(rotX_ * PId);
    view_.Z = temp.Z * cosX;
  }
  
  return view_;
}