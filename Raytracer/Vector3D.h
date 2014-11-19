#pragma once

#include <math.h>
#include "Matrix3x3D.h"

struct Vector3D
{
public:
	double X, Y, Z;

	Vector3D();
	Vector3D(double x, double y, double z);
	Vector3D(const double vec[]);

	inline double Dot(const Vector3D& other) { return Dot(*this, other); }
	inline Vector3D Cross(const Vector3D& other) { return Cross(*this, other); }
	inline double Length() { return sqrt(Dot(*this, *this)); }

	inline static double Dot(const Vector3D& lhs, const Vector3D& rhs) { return lhs.X * rhs.X + lhs.Y * rhs.Y + lhs.Z * rhs.Z; }
	inline static Vector3D Cross(const Vector3D& lhs, const Vector3D& rhs) { return Vector3D(lhs.Y * rhs.Z - lhs.Z * rhs.Y, lhs.Z * rhs.X - lhs.X * rhs.Z, lhs.X * rhs.Y - lhs.Y * rhs.X); }

	Vector3D& operator+=(const Vector3D& rhs);
	Vector3D& operator+=(const double& rhs);

	Vector3D& operator-=(const Vector3D& rhs);
	Vector3D& operator-=(const double& rhs);
	
	Vector3D& operator*=(const double& rhs);
	Vector3D& operator*=(const Matrix3x3D lhs);

	Vector3D& operator/=(const double& rhs);

private:
};

inline Vector3D operator+(Vector3D lhs, const Vector3D& rhs) { return lhs += rhs; }
inline Vector3D operator+(Vector3D lhs, const double& rhs) { return lhs += rhs; }
inline Vector3D operator+(const double& lhs, Vector3D rhs) { return rhs += lhs; }

inline Vector3D operator-(Vector3D lhs, const Vector3D& rhs) { return lhs -= rhs; }
inline Vector3D operator-(Vector3D lhs, const double& rhs) { return lhs -= rhs; }
inline Vector3D operator-(const double& lhs, Vector3D rhs) { return rhs -= lhs; }

inline Vector3D operator*(Vector3D lhs, const double& rhs) { return lhs *= rhs; }
inline Vector3D operator*(const double& lhs, Vector3D rhs) { return rhs *= lhs; }
inline Vector3D operator*(const Matrix3x3D& lhs, Vector3D rhs) { return rhs *= lhs; }

inline Vector3D operator/(Vector3D lhs, const double& rhs) { return lhs += rhs; }

inline bool operator==(const Vector3D& lhs, const Vector3D& rhs) { return lhs.X == rhs.X && lhs.Y == rhs.Y && rhs.Z == lhs.Z; }
inline bool operator!=(const Vector3D& lhs, const Vector3D& rhs) { return !operator==(lhs, rhs); }