#include <math.h>
#include <Matrix4x4D.h>
#include <Vector4D.h>

Matrix4x4D::Matrix4x4D()
{
  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++)
      this->elems[i][j] = 0;
}

Matrix4x4D::Matrix4x4D(const double elems[4][4])
{
  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++)
      this->elems[i][j] = elems[i][j];
}

Matrix4x4D::Matrix4x4D(const double elems[4])
{
  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++)
      this->elems[i][j] = elems[i * 4 + j];
}

Matrix4x4D::Matrix4x4D(const double& m00, const double& m01, const double& m02, const double& m03, const double& m10, const double& m11, const double& m12, const double& m13, const double& m20, const double& m21, const double& m22, const double& m23, const double& m30, const double& m31, const double& m32, const double& m33)
{
  elems[0][0] = m00;
  elems[0][1] = m01;
  elems[0][2] = m02;
  elems[0][3] = m03;
  elems[1][0] = m10;
  elems[1][1] = m11;
  elems[1][2] = m12;
  elems[1][3] = m13;
  elems[2][0] = m20;
  elems[2][1] = m21;
  elems[2][2] = m22;
  elems[2][3] = m23;
  elems[3][0] = m30;
  elems[3][1] = m31;
  elems[3][2] = m32;
  elems[3][3] = m33;
}

Matrix4x4D Matrix4x4D::Identity() { return Matrix4x4D(1, 0, 0, 0,
                                                      0, 1, 0, 0,
                                                      0, 0, 1, 0,
                                                      0, 0, 0, 1); }

void Matrix4x4D::Transpose()
{
  Matrix4x4D mat(*this);
  for (int i = 0; i < 4; ++i)
    for (int j = 0; j < 4; ++j)
      elems[i][j] = mat[j][i];
}

Matrix4x4D Matrix4x4D::Transposed() const
{
  Matrix4x4D mat(*this);
  for (int i = 0; i < 4; ++i)
    for (int j = 0; j < 4; ++j)
      mat[i][j] = elems[j][i];
  return mat;
}

Matrix4x4D Matrix4x4D::Transpose(const Matrix4x4D& matrix) { return matrix.Transposed(); }

double Matrix4x4D::Determinant()
{
  return elems[0][0] * elems[1][1] * elems[2][2] * elems[3][3] +
    elems[1][0] * elems[2][1] * elems[3][2] * elems[0][3] +
    elems[2][0] * elems[3][1] * elems[0][2] * elems[1][3] +
    elems[3][0] * elems[0][1] * elems[1][2] * elems[2][3] -
    elems[0][0] * elems[1][1] * elems[2][2] * elems[1][3] -
    elems[3][0] * elems[2][1] * elems[1][2] * elems[0][3] -
    elems[2][0] * elems[1][1] * elems[0][2] * elems[3][3] -
    elems[1][0] * elems[0][1] * elems[3][2] * elems[2][3];
}

void Matrix4x4D::Invert()
{
  double det = this->Determinant();
  if (!det)
    return;

  Matrix4x4D copy(*this);

  elems[0][0] = (copy[1][1] * copy[2][2] - copy[1][2] * copy[2][1]) / det;
  elems[0][1] = (copy[0][2] * copy[2][1] - copy[0][1] * copy[2][2]) / det;
  elems[0][2] = (copy[0][1] * copy[1][2] - copy[0][2] * copy[1][1]) / det;
  elems[1][0] = (copy[1][2] * copy[2][0] - copy[1][0] * copy[2][2]) / det;
  elems[1][1] = (copy[0][0] * copy[2][2] - copy[0][2] * copy[2][0]) / det;
  elems[1][2] = (copy[0][2] * copy[1][0] - copy[0][0] * copy[1][2]) / det;
  elems[2][0] = (copy[1][0] * copy[2][1] - copy[1][1] * copy[2][0]) / det;
  elems[2][1] = (copy[0][1] * copy[2][0] - copy[0][0] * copy[2][1]) / det;
  elems[2][2] = (copy[0][0] * copy[1][1] - copy[0][1] * copy[1][0]) / det;
}

Matrix4x4D Matrix4x4D::Inverse() const { Matrix4x4D m(*this); m.Invert(); return m; }
Matrix4x4D Matrix4x4D::Invert(const Matrix4x4D& matrix) { return matrix.Inverse(); }

double* Matrix4x4D::operator[](int i) { return elems[i]; }
const double* Matrix4x4D::operator[](int i) const { return elems[i]; }

Matrix4x4D& Matrix4x4D::operator*=(const Matrix4x4D& rhs)
{
  Matrix4x4D mat(*this);
  for (int i = 0; i < 4; ++i)
    for (int j = 0; j < 4; ++j)
      elems[i][j] = Vector4D::Dot(Vector4D(mat[i]), Vector4D(rhs[0][j], rhs[1][j], rhs[2][j], rhs[3][j]));

  return *this;
}

Matrix4x4D Matrix4x4D::createLookAt(Vector3D eye, Vector3D focus, Vector3D up)
{
  Vector3D forward = Vector3D::Normalise(eye - focus);
  Vector3D right = Vector3D::Normalise(Vector3D::Cross(up, forward));
  Vector3D newUp = Vector3D::Cross(forward, right);

  return Matrix4x4D(
    right.X,                      newUp.X,                      forward.X,                      0,
    right.Y,                      newUp.Y,                      forward.Y,                      0,
    right.Z,                      newUp.Z,                      forward.Z,                      0,
    -Vector3D::Dot(right, eye),   -Vector3D::Dot(newUp, eye),   -Vector3D::Dot(forward, eye),   1
  );
}

Matrix4x4D Matrix4x4D::createPerspectiveFieldOfView(double fieldOfView, double aspectRatio, double nearPlaneDistance, double farPlaneDistance)
{
  double yScale = 1.0 / tan(fieldOfView / 2);
  double xScale = yScale / aspectRatio;

  double planeDepth = nearPlaneDistance - farPlaneDistance;

  return Matrix4x4D(
    xScale,   0,        0,                                                  0,
    0,        yScale,   0,                                                  0,
    0,        0,        farPlaneDistance / planeDepth,                      -1,
    0,        0,        nearPlaneDistance * farPlaneDistance / planeDepth,  0
  );
}
