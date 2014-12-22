#pragma once

#include <Vector3D.h>

class Matrix4x4D
{
public:
  Matrix4x4D();
  Matrix4x4D(const double elems[4][4]);
  Matrix4x4D(const double elems[4]);
  Matrix4x4D(const double& m00, const double& m01, const double& m02, const double& m03,
    const double& m10, const double& m11, const double& m12, const double& m13,
    const double& m20, const double& m21, const double& m22, const double& m23,
    const double& m30, const double& m31, const double& m32, const double& m33);

  static Matrix4x4D Identity();

  void Transpose();
  Matrix4x4D Transposed() const;
  static Matrix4x4D Transpose(const Matrix4x4D& matrix);

  double Determinant();

  void Invert();
  Matrix4x4D Inverse() const;
  static Matrix4x4D Invert(const Matrix4x4D& matrix);

  double* operator[](int i);
  const double* operator[](int i) const;

  Matrix4x4D& operator*=(const Matrix4x4D& rhs);

  static Matrix4x4D createLookAt(Vector3D eye, Vector3D focus, Vector3D up);
  static Matrix4x4D createPerspectiveFieldOfView(double fieldOfView, double aspectRatio, double nearPlaneDistance, double farPlaneDistance);

private:
  double elems[4][4];
};

inline Matrix4x4D operator*(Matrix4x4D lhs, const Matrix4x4D& rhs) { return lhs *= rhs; }