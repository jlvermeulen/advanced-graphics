#include "Vector4D.h"

Vector4D::Vector4D() : X(0), Y(0), Z(0) { }
Vector4D::Vector4D(const double& x, const double& y, const double& z, const double& w) : X(x), Y(y), Z(z), W(w) { }
Vector4D::Vector4D(const double vec[]) { X = vec[0]; Y = vec[1]; Z = vec[2]; W = vec[3]; }

double Vector4D::Length() { return sqrt(Dot(*this, *this)); }

double Vector4D::Dot(const Vector4D& other) const { return Dot(*this, other); }
double Vector4D::Dot(const Vector4D& lhs, const Vector4D& rhs) { return lhs.X * rhs.X + lhs.Y * rhs.Y + lhs.Z * rhs.Z + lhs.W * rhs.W; }

void Vector4D::Normalise() { (*this) /= Length(); }
Vector4D Vector4D::Normalise(Vector4D vector) { vector.Normalise(); return vector; }

Vector4D& Vector4D::operator+=(const Vector4D& rhs)
{
  X += rhs.X;
  Y += rhs.Y;
  Z += rhs.Z;
  W += rhs.W;
  return *this;
}

Vector4D& Vector4D::operator+=(const double& rhs)
{
  X += rhs;
  Y += rhs;
  Z += rhs;
  W += rhs;
  return *this;
}

Vector4D& Vector4D::operator-=(const Vector4D& rhs)
{
  X -= rhs.X;
  Y -= rhs.Y;
  Z -= rhs.Z;
  W -= rhs.W;
  return *this;
}

Vector4D& Vector4D::operator-=(const double& rhs)
{
  X -= rhs;
  Y -= rhs;
  Z -= rhs;
  W -= rhs;
  return *this;
}

Vector4D& Vector4D::operator*=(const double& rhs)
{
  X *= rhs;
  Y *= rhs;
  Z *= rhs;
  W *= rhs;
  return *this;
}

Vector4D& Vector4D::operator*=(const Matrix4x4D& lhs)
{
  X = Dot(Vector4D(lhs[0]));
  Y = Dot(Vector4D(lhs[1]));
  Z = Dot(Vector4D(lhs[2]));
  W = Dot(Vector4D(lhs[3]));
  return *this;
}

Vector4D& Vector4D::operator/=(const double& rhs)
{
  X /= rhs;
  Y /= rhs;
  Z /= rhs;
  W /= rhs;
  return *this;
}

double& Vector4D::operator[](int index) { return index == 0 ? X : index == 1 ? Y : index == 2 ? Z : W; }
const double& Vector4D::operator[](int index) const { return index == 0 ? X : index == 1 ? Y : index == 2 ? Z : W; };