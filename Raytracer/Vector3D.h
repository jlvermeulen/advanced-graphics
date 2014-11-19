#pragma once

#include <math.h>
#include "Matrix3x3D.h"

struct Vector3D
{
public:
	double X, Y, Z;
	Vector3D() { X = Y = Z = 0; }
	Vector3D(double x, double y, double z) : X(x), Y(y), Z(z) { }
	Vector3D(const double vec[]) { X = vec[0]; Y = vec[1]; Z = vec[2]; }
	Vector3D(const Vector3D& vec) { X = vec.X; Y = vec.Y; Z = vec.Z; }
	~Vector3D() { }

	inline double Dot(const Vector3D& other) { return Dot(*this, other); }
	inline Vector3D Cross(const Vector3D& other) { return Cross(*this, other); }
	inline double Length() { return sqrt(Dot(*this, *this)); }

	inline static double Dot(const Vector3D& lhs, const Vector3D& rhs) { return lhs.X * rhs.X + lhs.Y * rhs.Y + lhs.Z * rhs.Z; }
	inline static Vector3D Cross(const Vector3D& lhs, const Vector3D& rhs) { return Vector3D(lhs.Y * rhs.Z - lhs.Z * rhs.Y, lhs.Z * rhs.X - lhs.X * rhs.Z, lhs.X * rhs.Y - lhs.Y * rhs.X); }

	Vector3D& operator+=(const Vector3D& rhs)
	{
		X += rhs.X;
		Y += rhs.Y;
		Z += rhs.Z;
		return *this;
	}

	Vector3D& operator+=(const double& rhs)
	{
		X += rhs;
		Y += rhs;
		Z += rhs;
		return *this;
	}

	Vector3D& operator-=(const Vector3D& rhs)
	{
		X -= rhs.X;
		Y -= rhs.Y;
		Z -= rhs.Z;
		return *this;
	}

	Vector3D& operator-=(const double& rhs)
	{
		X -= rhs;
		Y -= rhs;
		Z -= rhs;
		return *this;
	}
	
	Vector3D& operator*=(const double& rhs)
	{
		X *= rhs;
		Y *= rhs;
		Z *= rhs;
		return *this;
	}

	Vector3D& operator*=(const Matrix3x3D lhs)
	{
		X = Dot(Vector3D(lhs[0]));
		Y = Dot(Vector3D(lhs[1]));
		Z = Dot(Vector3D(lhs[2]));
		return *this;
	}

	Vector3D& operator/=(const double& rhs)
	{
		X /= rhs;
		Y /= rhs;
		Z /= rhs;
		return *this;
	}

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