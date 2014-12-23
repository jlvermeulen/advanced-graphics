#pragma once

#include <Vector3D.h>

class Camera
{
public:
  Camera();
  ~Camera();

public:
  void MoveForward(const float& distance);
  void MoveRight(const float& distance);
  void MoveUpward(const float& distance);
  void RotateX(float angle);
  void RotateY(float angle);
  void RotateZ(float angle);


  Vector3D getEye() const { return eye_; }
  Vector3D getFocus() const { return focus_; }
  Vector3D getUp() const { return up_;  }

private:
  Vector3D eye_;
  Vector3D focus_;
  Vector3D right_;
  Vector3D up_;

  float rotX_;
  float rotY_;
  float rotZ_;
};