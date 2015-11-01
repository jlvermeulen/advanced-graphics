#pragma once

#include "Vector3F.h"

#include <math.h>

class Quaternion
{
public:
  float X, Y, Z, W;

  Quaternion();
  Quaternion(Vector3F xyz, float w);
  Quaternion(float x, float y, float z, float w);

  float length() const;

  void normalize();
  static Quaternion normalize(Quaternion quat);

  void conjugate();
  static Quaternion conjugate(Quaternion quat);

  Quaternion& operator*=(const Quaternion& rhs);
  Quaternion& operator/=(const float& rhs);
};

inline Quaternion operator*(Quaternion lhs, const Quaternion& rhs) { return lhs *= rhs; }
inline Quaternion operator/(Quaternion lhs, const float& rhs) { return lhs /= rhs; }