#pragma once

#include "Vector3F.h"
#include "Quaternion.h"

class Camera
{
public:
  Camera();
  Camera(Vector3F eye, Vector3F focus, Vector3F up, int width = 1280, int height = 720, float fovY = 40.0f, float zNear = 0.1f, float zFar = 1000.0f, float aperture = 3.0f, float focalDepth = 4.0f);
  ~Camera();

public:
  void MoveForward(const float& distance);
  void MoveRight(const float& distance);
  void MoveUpward(const float& distance);
  void Rotate(const Vector3F& axis, float angle);
  void RotateX(float angle);
  void RotateY(float angle);
  void RotateZ(float angle);

  void setX(float x) { eye_.X = x; }
  void setY(float y) { eye_.Y = y; }
  void setZ(float z) { eye_.Z = z; }

  Vector3F Eye() const { return eye_; }
  Vector3F Focus() const { return focus_; }
  Vector3F Right() const { return right_; }
  Vector3F Up() const { return up_;  }
  float FovY() const { return fovY_; }
  float ZNear() const { return zNear_; }
  float ZFar() const { return zFar_; }

public:
  int Width;
  int Height;

  float Aperture;
  float FocalDepth;

private:
  Vector3F eye_;
  Vector3F focus_;
  Vector3F right_;
  Vector3F up_;

  Quaternion orientation_;
  
  float fovY_;
  float zNear_;
  float zFar_;
};