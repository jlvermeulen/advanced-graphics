#pragma once

#include <math.h>
#include "Matrix3x3D.h"

class Vector3D
{
public:
	double X, Y, Z;

	Vector3D();
	Vector3D(const double& x, const double& y, const double& z);
	Vector3D(const double vec[]);

	double Length() const;
	double LengthSquared() const;

	double Dot(const Vector3D& other) const;
	static double Dot(const Vector3D& lhs, const Vector3D& rhs);

	Vector3D Cross(const Vector3D& other) const;
	static Vector3D Cross(const Vector3D& lhs, const Vector3D& rhs);

	void Normalise();
	static Vector3D Normalise(Vector3D vector);

	void Reflect(const Vector3D& normal);
	Vector3D Reflected(const Vector3D& normal) const;
	static Vector3D Reflect(Vector3D incident, const Vector3D& normal);

	void Refract(const Vector3D& normal, const double& n1, const double& n2);
	Vector3D Refracted(const Vector3D& normal, const double& n1, const double& n2) const;
	static Vector3D Refract(Vector3D incident, const Vector3D& normal, const double& n1, const double& n2);

	Vector3D& operator+=(const Vector3D& rhs);
	Vector3D& operator+=(const double& rhs);

	Vector3D& operator-=(const Vector3D& rhs);
	Vector3D& operator-=(const double& rhs);
	
	Vector3D& operator*=(const double& rhs);
	Vector3D& operator*=(const Matrix3x3D& lhs);

	Vector3D& operator/=(const double& rhs);

	double& operator[](int index);
	const double& operator[](int index) const;

public:
	static Vector3D Forward;
	static Vector3D Right;
	static Vector3D Up;

private:
};

inline Vector3D operator+(Vector3D lhs, const Vector3D& rhs) { return lhs += rhs; }
inline Vector3D operator+(Vector3D lhs, const double& rhs) { return lhs += rhs; }
inline Vector3D operator+(const double& lhs, Vector3D rhs) { return rhs += lhs; }

inline Vector3D operator-(Vector3D lhs, const Vector3D& rhs) { return lhs -= rhs; }
inline Vector3D operator-(Vector3D lhs, const double& rhs) { return lhs -= rhs; }
inline Vector3D operator-(const double& lhs, Vector3D rhs) { return rhs -= lhs; }
inline Vector3D operator-(Vector3D rhs) { return rhs *= -1; }

inline Vector3D operator*(Vector3D lhs, const double& rhs) { return lhs *= rhs; }
inline Vector3D operator*(const double& lhs, Vector3D rhs) { return rhs *= lhs; }
inline Vector3D operator*(const Matrix3x3D& lhs, Vector3D rhs) { return rhs *= lhs; }

inline Vector3D operator/(Vector3D lhs, const double& rhs) { return lhs /= rhs; }

inline bool operator==(const Vector3D& lhs, const Vector3D& rhs) { return lhs.X == rhs.X && lhs.Y == rhs.Y && lhs.Z == rhs.Z; }
inline bool operator!=(const Vector3D& lhs, const Vector3D& rhs) { return !operator==(lhs, rhs); }