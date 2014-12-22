#pragma once

#include <Matrix4x4D.h>
#include <Vector3D.h>

class Camera
{
public:
  Camera(Vector3D eye, Vector3D focus, Vector3D up, double aspectRatio = 16.0f / 9.0f);
  ~Camera();

private:
  void updateViewMatrix();

  Vector3D getEye() { return eye_; }
  void setEye(Vector3D eye);

  Vector3D getFocus() { return focus_; }
  void setFocus(Vector3D focus);

  Vector3D getUp() { return up_; }
  void setUp(Vector3D up);

private:
  Matrix4x4D view_;
  Matrix4x4D projection_;

  Vector3D eye_;
  Vector3D focus_;
  Vector3D up_;
};