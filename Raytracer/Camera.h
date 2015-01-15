#pragma once

#include <Vector3D.h>

class Camera
{
public:
  Camera();
  Camera(Vector3D eye, Vector3D focus, Vector3D up, int width = 1280, int height = 720, float fovY = 40.0, float zNear = 0.1, float zFar = 1000.0);
  ~Camera();

public:
  void MoveForward(const float& distance);
  void MoveRight(const float& distance);
  void MoveUpward(const float& distance);
  void RotateX(float angle);
  void RotateY(float angle);
  void RotateZ(float angle);

  Vector3D Eye() const { return eye_; }
  Vector3D Focus() const { return focus_; }
  Vector3D Right() const { return right_; }
  Vector3D Up() const { return up_;  }
  float FovY() const { return fovY_; }
  float ZNear() const { return zNear_; }
  float ZFar() const { return zFar_; }

public:
  int Width;
  int Height;

private:
  Vector3D eye_;
  Vector3D focus_;
  Vector3D right_;
  Vector3D up_;

  float rotX_;
  float rotY_;
  float rotZ_;
  
  float fovY_;
  float zNear_;
  float zFar_;
};