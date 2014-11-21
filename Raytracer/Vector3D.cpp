#include "Vector3D.h"

Vector3D::Vector3D(const double& x, const double& y, const double& z) : X(x), Y(y), Z(z) { }
Vector3D::Vector3D(const double vec[]) { X = vec[0]; Y = vec[1]; Z = vec[2]; }

double Vector3D::Length() { return sqrt(Dot(*this, *this)); }

double Vector3D::Dot(const Vector3D& other) const { return Dot(*this, other); }
double Vector3D::Dot(const Vector3D& lhs, const Vector3D& rhs) { return lhs.X * rhs.X + lhs.Y * rhs.Y + lhs.Z * rhs.Z; }

Vector3D Vector3D::Cross(const Vector3D& other) const { return Cross(*this, other); }
Vector3D Vector3D::Cross(const Vector3D& lhs, const Vector3D& rhs) { return Vector3D(lhs.Y * rhs.Z - lhs.Z * rhs.Y, lhs.Z * rhs.X - lhs.X * rhs.Z, lhs.X * rhs.Y - lhs.Y * rhs.X); }

void Vector3D::Normalise() { (*this) /= Length(); }
Vector3D Vector3D::Normalise(Vector3D vector) { vector.Normalise(); return vector; }

void Vector3D::Reflect(const Vector3D& normal) { (*this) += 2 * Dot(normal, -(*this)) * normal; }
Vector3D Vector3D::Reflected(const Vector3D& normal) const { return Reflect((*this), normal); }
Vector3D Vector3D::Reflect(Vector3D incident, const Vector3D& normal) { incident.Reflect(normal); return incident; }

void Vector3D::Refract(const Vector3D& normal, const double& nFrom, const double& nTo)
{
	double r = nFrom / nTo;
	double c = -Dot((*this), normal);

	double radicand = 1 - r * r * (1 - c * c);
	if (radicand < 0)
		Reflect(normal);

	(*this) *= r;
	(*this) += (r * c - sqrt(radicand)) * normal;
}
Vector3D Vector3D::Refracted(const Vector3D& normal, const double& nFrom, const double& nTo) const { return Refract((*this), normal, nFrom, nTo); }
Vector3D Vector3D::Refract(Vector3D incident, const Vector3D& normal, const double& nFrom, const double& nTo) { incident.Refract(normal, nFrom, nTo); return incident; }

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