#include "Vector3D.h"

//---------------------------------------------------------------------------------------
Vector3D::Vector3D() :
	X(0),
	Y(0),
	Z(0)
{ }

//---------------------------------------------------------------------------------------
Vector3D::Vector3D(const double& x, const double& y, const double& z) :
	X(x),
	Y(y),
	Z(z)
{ }

//---------------------------------------------------------------------------------------
Vector3D::Vector3D(const double vec[]) :
	X(vec[0]),
	Y(vec[1]),
	Z(vec[2])
{ }

//---------------------------------------------------------------------------------------
Vector3D::Vector3D(const Vector3D& rhs) :
	X(rhs.X),
	Y(rhs.Y),
	Z(rhs.Z)
{ }

//---------------------------------------------------------------------------------------
Vector3D::~Vector3D()
{ }

//---------------------------------------------------------------------------------------
double Vector3D::Length() const
{
	return sqrt(X * X + Y * Y + Z * Z);
}

//---------------------------------------------------------------------------------------
double Vector3D::LengthSquared() const
{
	return (X * X + Y * Y + Z * Z);
}

//---------------------------------------------------------------------------------------
double Vector3D::Distance(const Vector3D& v) const
{
	double dx = v.X - X;
	double dy = v.Y - Y;
	double dz = v.Z - Z;

	return sqrt(dx * dx + dy * dy + dz * dz);
}

//---------------------------------------------------------------------------------------
double Vector3D::DistanceSquared(const Vector3D& v) const
{
	double dx = v.X - X;
	double dy = v.Y - Y;
	double dz = v.Z - Z;

	return dx * dx + dy * dy + dz * dz;
}

//---------------------------------------------------------------------------------------
double Vector3D::Dot(const Vector3D& v) const
{
	return (X * v.X + Y * v.Y + Z * v.Z);
}

//---------------------------------------------------------------------------------------
double Vector3D::Dot(const Vector3D& v1, const Vector3D& v2)
{
	return (v1.X * v2.X + v1.Y * v2.Y + v1.Z * v2.Z);
}

//---------------------------------------------------------------------------------------
Vector3D Vector3D::Cross(const Vector3D& v) const
{
	return Vector3D(
		Y * v.Z - Z * v.Y,
		Z * v.X - X * v.Z,
		X * v.Y - Y * v.X
	);
}

//---------------------------------------------------------------------------------------
Vector3D Vector3D::Cross(const Vector3D& v1, const Vector3D& v2) {
	return Vector3D(
		v1.Y * v2.Z - v1.Z * v2.Y,
		v1.Z * v2.X - v1.X * v2.Z,
		v1.X * v2.Y - v1.Y * v2.X
	);
}

//---------------------------------------------------------------------------------------
void Vector3D::Normalise()
{
	double n = LengthSquared();

	// Not normalized yet.
	if (n != 1.0)
	{
		n = sqrt(n);

		*this /= n;
	}
}

//---------------------------------------------------------------------------------------
Vector3D Vector3D::Normalise(Vector3D v) {
	Vector3D result(v);

	result.Normalise();

	return result;
}