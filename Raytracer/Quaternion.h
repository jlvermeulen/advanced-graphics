#pragma once

#include "Vector3D.h"

#include <math.h>

class Quaternion
{
public:
  double X, Y, Z, W;

  Quaternion();
  Quaternion(Vector3D xyz, double w);
  Quaternion(double x, double y, double z, double w);

  double length() const;

  void normalize();
  static Quaternion normalize(Quaternion quat);

  void conjugate();
  static Quaternion conjugate(Quaternion quat);

  Quaternion& operator*=(const Quaternion& rhs);
  Quaternion& operator/=(const double& rhs);
};

inline Quaternion operator*(Quaternion lhs, const Quaternion& rhs) { return lhs *= rhs; }
inline Quaternion operator/(Quaternion lhs, const double& rhs) { return lhs /= rhs; }