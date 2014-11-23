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

	inline Vector3D operator+(const Vector3D& rhs);
	inline Vector3D operator+(const double& rhs);
	friend inline Vector3D operator+(const double& lhs, const Vector3D& rhs);

	inline Vector3D operator-();
	inline Vector3D operator-(const Vector3D& rhs);
	inline Vector3D operator-(const double& rhs);
	friend inline Vector3D operator-(const double& lhs, const Vector3D& rhs);

	inline Vector3D operator*(const double& rhs);
	friend inline Vector3D operator*(const double& lhs, const Vector3D& rhs);
	friend inline Vector3D operator*(const Matrix3x3D& lhs, const Vector3D& rhs);
	
	inline Vector3D operator/(const double& rhs);
	
	inline void operator+=(const Vector3D& rhs);
	inline void operator+=(const double& rhs);

	inline void operator-=(const Vector3D& rhs);
	inline void operator-=(const double& rhs);

	inline void operator*=(const double& rhs);
	inline void operator/=(const double& rhs);

	inline bool operator==(const Vector3D& rhs);
	inline bool operator!=(const Vector3D& rhs);

private:
};

#include "Vector3D.inl"