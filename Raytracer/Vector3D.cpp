#include "Vector3D.h"

Vector3D::Vector3D() : X(0), Y(0), Z(0) { }
Vector3D::Vector3D(const double& x, const double& y, const double& z) : X(x), Y(y), Z(z) { }
Vector3D::Vector3D(const double vec[]) { X = vec[0]; Y = vec[1]; Z = vec[2]; }

double Vector3D::Length() { return sqrt(LengthSquared()); }
double Vector3D::LengthSquared() { return Dot(*this, *this); }

double Vector3D::Dot(const Vector3D& other) const { return Dot(*this, other); }
double Vector3D::Dot(const Vector3D& lhs, const Vector3D& rhs) { return lhs.X * rhs.X + lhs.Y * rhs.Y + lhs.Z * rhs.Z; }

Vector3D Vector3D::Cross(const Vector3D& other) const { return Cross(*this, other); }
Vector3D Vector3D::Cross(const Vector3D& lhs, const Vector3D& rhs) { return Vector3D(lhs.Y * rhs.Z - lhs.Z * rhs.Y, lhs.Z * rhs.X - lhs.X * rhs.Z, lhs.X * rhs.Y - lhs.Y * rhs.X); }

void Vector3D::Normalise() { (*this) /= Length(); }
Vector3D Vector3D::Normalise(Vector3D vector) { vector.Normalise(); return vector; }

Vector3D& Vector3D::operator+=(const Vector3D& rhs)
{
	X += rhs.X;
	Y += rhs.Y;
	Z += rhs.Z;
	return *this;
}

Vector3D& Vector3D::operator+=(const double& rhs)
{
	X += rhs;
	Y += rhs;
	Z += rhs;
	return *this;
}

Vector3D& Vector3D::operator-=(const Vector3D& rhs)
{
	X -= rhs.X;
	Y -= rhs.Y;
	Z -= rhs.Z;
	return *this;
}

Vector3D& Vector3D::operator-=(const double& rhs)
{
	X -= rhs;
	Y -= rhs;
	Z -= rhs;
	return *this;
}

Vector3D& Vector3D::operator*=(const double& rhs)
{
	X *= rhs;
	Y *= rhs;
	Z *= rhs;
	return *this;
}

Vector3D& Vector3D::operator*=(const Matrix3x3D& lhs)
{
	X = Dot(Vector3D(lhs[0]));
	Y = Dot(Vector3D(lhs[1]));
	Z = Dot(Vector3D(lhs[2]));
	return *this;
}

Vector3D& Vector3D::operator/=(const double& rhs)
{
	X /= rhs;
	Y /= rhs;
	Z /= rhs;
	return *this;
}

double& Vector3D::operator[](int index) { return index == 0 ? X : index == 1 ? Y : Z; }
const double& Vector3D::operator[](int index) const { return index == 0 ? X : index == 1 ? Y : Z; };

Vector3D Vector3D::Forward = Vector3D(0, 0, -1);
Vector3D Vector3D::Right = Vector3D(1, 0, 0);
Vector3D Vector3D::Up = Vector3D(0, 1, 0);