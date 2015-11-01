#include "Quaternion.h"

//--------------------------------------------------------------------------------
Quaternion::Quaternion() :
  X(0),
  Y(0),
  Z(0),
  W(0)
{
}
//--------------------------------------------------------------------------------
Quaternion::Quaternion(Vector3F xyz, float w) :
  X(xyz.X),
  Y(xyz.Y),
  Z(xyz.Z),
  W(w)
{

}

//--------------------------------------------------------------------------------
Quaternion::Quaternion(float x, float y, float z, float w) :
  X(x),
  Y(y),
  Z(z),
  W(w)
{

}

//--------------------------------------------------------------------------------
float Quaternion::length() const { return sqrtf(X * X + Y * Y + Z * Z + W * W); }

//--------------------------------------------------------------------------------
void Quaternion::normalize()
{
  *this /= length();
}

//--------------------------------------------------------------------------------
Quaternion Quaternion::normalize(Quaternion quat)
{
  quat.normalize();

  return quat;
}

//--------------------------------------------------------------------------------
void Quaternion::conjugate()
{
  X *= -1;
  Y *= -1;
  Z *= -1;
}

//--------------------------------------------------------------------------------
Quaternion Quaternion::conjugate(Quaternion quat)
{
  quat.conjugate();

  return quat;
}

//--------------------------------------------------------------------------------
Quaternion& Quaternion::operator*=(const Quaternion& rhs)
{
  Quaternion lhs(*this);

  X = lhs.W * rhs.X + lhs.X * rhs.W + lhs.Y * rhs.Z - lhs.Z * rhs.Y;
  Y = lhs.W * rhs.Y - lhs.X * rhs.Z + lhs.Y * rhs.W + lhs.Z * rhs.X;
  Z = lhs.W * rhs.Z + lhs.X * rhs.Y - lhs.Y * rhs.X + lhs.Z * rhs.W;
  W = lhs.W * rhs.W - lhs.X * rhs.X - lhs.Y * rhs.Y - lhs.Z * rhs.Z;

  return *this;
}

//--------------------------------------------------------------------------------
Quaternion& Quaternion::operator/=(const float& rhs)
{
  X /= rhs;
  Y /= rhs;
  Z /= rhs;
  W /= rhs;

  return *this;
}