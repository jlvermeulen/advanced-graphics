#pragma once

#include <math.h>
#include "Matrix3x3D.h"

struct Vector3D
{
public:
	double X, Y, Z;

	Vector3D();
	Vector3D(const double& x, const double& y, const double& z);
	Vector3D(const double vec[]);
	Vector3D(const Vector3D& rhs);

	~Vector3D();

	double Length();

	double Dot(const Vector3D& other) const;
	static double Dot(const Vector3D& lhs, const Vector3D& rhs);

	Vector3D Cross(const Vector3D& other) const;
	static Vector3D Cross(const Vector3D& lhs, const Vector3D& rhs);

	void Normalise();
	static Vector3D Normalise(Vector3D vector);

	void Reflect(const Vector3D& normal);
	Vector3D Reflected(const Vector3D& normal) const;
	static Vector3D Reflect(Vector3D incident, const Vector3D& normal);

	void Refract(const Vector3D& normal, const double& nFrom, const double& nTo);
	Vector3D Refracted(const Vector3D& normal, const double& nFrom, const double& nTo) const;
	static Vector3D Refract(Vector3D incident, const Vector3D& normal, const double& nFrom, const double& nTo);

	Vector3D& operator+=(const Vector3D& rhs);
	Vector3D& operator+=(const double& rhs);

	Vector3D& operator-=(const Vector3D& rhs);
	Vector3D& operator-=(const double& rhs);
	
	Vector3D& operator*=(const double& rhs);
	Vector3D& operator*=(const Matrix3x3D lhs);

	Vector3D& operator/=(const double& rhs);

private:
};

inline Vector3D operator+(const Vector3D& lhs, const Vector3D& rhs)
{
	Vector3D result(lhs);

	result.X += rhs.X;
	result.Y += rhs.Y;
	result.Z += rhs.Z;

	return result;
}

inline Vector3D operator+(const Vector3D& lhs, const double& rhs)
{
	Vector3D result(lhs);

	result.X += rhs;
	result.Y += rhs;
	result.Z += rhs;

	return result;
}

inline Vector3D operator+(const double& lhs, const Vector3D& rhs)
{
	Vector3D result(rhs);

	result.X += lhs;
	result.Y += lhs;
	result.Z += lhs;

	return result;
}

inline Vector3D operator-(const Vector3D& lhs, const Vector3D& rhs)
{
	Vector3D result(lhs);

	result.X -= rhs.X;
	result.Y -= rhs.Y;
	result.Z -= rhs.Z;

	return result;
}

inline Vector3D operator-(const Vector3D& lhs, const double& rhs)
{
	Vector3D result(lhs);

	result.X -= rhs;
	result.Y -= rhs;
	result.Z -= rhs;

	return result;
}

inline Vector3D operator-(const double& lhs, const Vector3D& rhs)
{
	Vector3D result(rhs);

	result.X -= lhs;
	result.Y -= lhs;
	result.Z -= lhs;

	return result;
}

inline Vector3D operator-(const Vector3D& rhs)
{
	Vector3D result(rhs);

	result.X *= -1;
	result.Y *= -1;
	result.Z *= -1;

	return result;
}

inline Vector3D operator*(const Vector3D& lhs, const double& rhs)
{
	Vector3D result(lhs);

	result.X *= rhs;
	result.Y *= rhs;
	result.Z *= rhs;

	return result;
}

inline Vector3D operator*(const double& lhs, const Vector3D& rhs)
{
	Vector3D result(rhs);

	result.X *= lhs;
	result.Y *= lhs;
	result.Z *= lhs;

	return result;
}

inline Vector3D operator*(const Matrix3x3D& lhs, const Vector3D& rhs)
{
	Vector3D result(rhs);

	result.X = result.Dot(Vector3D(lhs[0]));
	result.Y = result.Dot(Vector3D(lhs[1]));
	result.Z = result.Dot(Vector3D(lhs[2]));

	return result;
}

inline Vector3D operator/(const Vector3D& lhs, const double& rhs)
{
	Vector3D result(lhs);

	result.X /= rhs;
	result.Y /= rhs;
	result.Z /= rhs;

	return result;
}

inline bool operator==(const Vector3D& lhs, const Vector3D& rhs)
{
	return lhs.X == rhs.X && lhs.Y == rhs.Y && lhs.Z == rhs.Z;
}

inline bool operator!=(const Vector3D& lhs, const Vector3D& rhs)
{
	return lhs.X != rhs.X || lhs.Y != rhs.Y || lhs.Z != rhs.Z;
}