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

	double Length() const;
	double LengthSquared() const;

	double Distance(const Vector3D& v) const;
	double DistanceSquared(const Vector3D& v) const;

	double Dot(const Vector3D& v) const;
	static double Dot(const Vector3D& v1, const Vector3D& v2);

	Vector3D Cross(const Vector3D& other) const;
	static Vector3D Cross(const Vector3D& v1, const Vector3D& v2);

	void Normalise();
	static Vector3D Normalise(Vector3D vector);

	inline Vector3D operator+(const Vector3D& rhs) const;
	inline Vector3D operator+(const double& rhs) const;
	friend inline Vector3D operator+(const double& lhs, const Vector3D& rhs);

	inline Vector3D operator-() const;
	inline Vector3D operator-(const Vector3D& rhs) const;
	inline Vector3D operator-(const double& rhs) const;
	friend inline Vector3D operator-(const double& lhs, const Vector3D& rhs);

	inline Vector3D operator*(const double& rhs) const;
	friend inline Vector3D operator*(const double& lhs, const Vector3D& rhs);
	friend inline Vector3D operator*(const Matrix3x3D& lhs, const Vector3D& rhs);
	
	inline Vector3D operator/(const double& rhs) const;
	
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