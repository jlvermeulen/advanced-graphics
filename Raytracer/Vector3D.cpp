#pragma once

#include "Vector3D.h"

Vector3D::Vector3D() : X(0), Y(0), Z(0) { }
Vector3D::Vector3D(double x, double y, double z) : X(x), Y(y), Z(z) { }
Vector3D::Vector3D(const double vec[]) { X = vec[0]; Y = vec[1]; Z = vec[2]; }

Vector3D& Vector3D::operator+=(const Vector3D& rhs)
{
	Vector3D v = Vector3D() * 5;
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

Vector3D& Vector3D::operator*=(const Matrix3x3D lhs)
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