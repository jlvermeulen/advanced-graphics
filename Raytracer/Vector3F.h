#pragma once

#include <math.h>
#include "Matrix3x3F.h"

class Vector3F
{
public:
	float X, Y, Z;

	Vector3F();
	Vector3F(const float& x, const float& y, const float& z);
	Vector3F(const float vec[]);

	float Length() const;
	float LengthSquared() const;

	float Dot(const Vector3F& other) const;
	static float Dot(const Vector3F& lhs, const Vector3F& rhs);

	Vector3F Cross(const Vector3F& other) const;
	static Vector3F Cross(const Vector3F& lhs, const Vector3F& rhs);

	void Normalise();
	static Vector3F Normalise(Vector3F vector);

	void Reflect(const Vector3F& normal);
	Vector3F Reflected(const Vector3F& normal) const;
	static Vector3F Reflect(Vector3F incident, const Vector3F& normal);

	void Refract(const Vector3F& normal, const float& n1, const float& n2);
	Vector3F Refracted(const Vector3F& normal, const float& n1, const float& n2) const;
	static Vector3F Refract(Vector3F incident, const Vector3F& normal, const float& n1, const float& n2);

	Vector3F& operator+=(const Vector3F& rhs);
	Vector3F& operator+=(const float& rhs);

	Vector3F& operator-=(const Vector3F& rhs);
	Vector3F& operator-=(const float& rhs);
	
	Vector3F& operator*=(const float& rhs);
	Vector3F& operator*=(const Matrix3x3F& lhs);

	Vector3F& operator/=(const float& rhs);

	float& operator[](int index);
	const float& operator[](int index) const;

public:
	static Vector3F Forward;
	static Vector3F Right;
	static Vector3F Up;

private:
};

inline Vector3F operator+(Vector3F lhs, const Vector3F& rhs) { return lhs += rhs; }
inline Vector3F operator+(Vector3F lhs, const float& rhs) { return lhs += rhs; }
inline Vector3F operator+(const float& lhs, Vector3F rhs) { return rhs += lhs; }

inline Vector3F operator-(Vector3F lhs, const Vector3F& rhs) { return lhs -= rhs; }
inline Vector3F operator-(Vector3F lhs, const float& rhs) { return lhs -= rhs; }
inline Vector3F operator-(const float& lhs, Vector3F rhs) { return rhs -= lhs; }
inline Vector3F operator-(Vector3F rhs) { return rhs *= -1; }

inline Vector3F operator*(Vector3F lhs, const float& rhs) { return lhs *= rhs; }
inline Vector3F operator*(const float& lhs, Vector3F rhs) { return rhs *= lhs; }
inline Vector3F operator*(const Matrix3x3F& lhs, Vector3F rhs) { return rhs *= lhs; }

inline Vector3F operator/(Vector3F lhs, const float& rhs) { return lhs /= rhs; }

inline bool operator==(const Vector3F& lhs, const Vector3F& rhs) { return lhs.X == rhs.X && lhs.Y == rhs.Y && lhs.Z == rhs.Z; }
inline bool operator!=(const Vector3F& lhs, const Vector3F& rhs) { return !operator==(lhs, rhs); }
