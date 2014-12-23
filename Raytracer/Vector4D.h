#pragma once

#include <math.h>
#include <Matrix4x4D.h>

class Vector4D
{
public:
  double X, Y, Z, W;

  Vector4D();
  Vector4D(const double& x, const double& y, const double& z, const double& w);
  Vector4D(const double vec[]);

  double Length();

  double Dot(const Vector4D& other) const;
  static double Dot(const Vector4D& lhs, const Vector4D& rhs);

  void Normalise();
  static Vector4D Normalise(Vector4D vector);

  Vector4D& operator+=(const Vector4D& rhs);
  Vector4D& operator+=(const double& rhs);

  Vector4D& operator-=(const Vector4D& rhs);
  Vector4D& operator-=(const double& rhs);

  Vector4D& operator*=(const double& rhs);
  Vector4D& operator*=(const Matrix4x4D& lhs);

  Vector4D& operator/=(const double& rhs);

  double& operator[](int index);
  const double& operator[](int index) const;

private:
};

inline Vector4D operator+(Vector4D lhs, const Vector4D& rhs) { return lhs += rhs; }
inline Vector4D operator+(Vector4D lhs, const double& rhs) { return lhs += rhs; }
inline Vector4D operator+(const double& lhs, Vector4D rhs) { return rhs += lhs; }

inline Vector4D operator-(Vector4D lhs, const Vector4D& rhs) { return lhs -= rhs; }
inline Vector4D operator-(Vector4D lhs, const double& rhs) { return lhs -= rhs; }
inline Vector4D operator-(const double& lhs, Vector4D rhs) { return rhs -= lhs; }
inline Vector4D operator-(Vector4D rhs) { return rhs *= -1; }

inline Vector4D operator*(Vector4D lhs, const double& rhs) { return lhs *= rhs; }
inline Vector4D operator*(const double& lhs, Vector4D rhs) { return rhs *= lhs; }
inline Vector4D operator*(const Matrix4x4D& lhs, Vector4D rhs) { return rhs *= lhs; }

inline Vector4D operator/(Vector4D lhs, const double& rhs) { return lhs /= rhs; }

inline bool operator==(const Vector4D& lhs, const Vector4D& rhs) { return lhs.X == rhs.X && lhs.Y == rhs.Y && rhs.Z == lhs.Z && lhs.W == rhs.W; }
inline bool operator!=(const Vector4D& lhs, const Vector4D& rhs) { return !operator==(lhs, rhs); }