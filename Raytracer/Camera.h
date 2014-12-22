#pragma once

#include <Vector3D.h>

class Camera
{
public:
  Camera();
  ~Camera();

public:
  void MoveForward(const float& distance);
  void MoveSideways(const float& distance);
  void Move(const Vector3D& direction);
  void RotateX(float angle);
  void RotateY(float angle);
  void RotateZ(float angle);

  Vector3D getPosition() const { return eye_; }
  Vector3D getRotation() const { return Vector3D(rotX_, rotY_, rotZ_); }
  const Vector3D& getViewDirection();

private:
  bool changed_;

  Vector3D eye_;
  Vector3D focus_;
  Vector3D view_;

  float rotX_;
  float rotY_;
  float rotZ_;
};