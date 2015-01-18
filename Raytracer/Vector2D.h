#pragma once

#include <math.h>

class Vector2D
{
public:
	double X, Y;

	Vector2D();
	Vector2D(const double& x, const double& y);
	Vector2D(const double vec[]);

	double Length();
  double LengthSquared();

	double Dot(const Vector2D& other) const;
	static double Dot(const Vector2D& lhs, const Vector2D& rhs);

	void Normalise();
	static Vector2D Normalise(Vector2D vector);

	Vector2D& operator+=(const Vector2D& rhs);
	Vector2D& operator+=(const double& rhs);

	Vector2D& operator-=(const Vector2D& rhs);
	Vector2D& operator-=(const double& rhs);
	
	Vector2D& operator*=(const double& rhs);

	Vector2D& operator/=(const double& rhs);

	double& operator[](int index);
	const double& operator[](int index) const;
};

inline Vector2D operator+(Vector2D lhs, const Vector2D& rhs) { return lhs += rhs; }
inline Vector2D operator+(Vector2D lhs, const double& rhs) { return lhs += rhs; }
inline Vector2D operator+(const double& lhs, Vector2D rhs) { return rhs += lhs; }

inline Vector2D operator-(Vector2D lhs, const Vector2D& rhs) { return lhs -= rhs; }
inline Vector2D operator-(Vector2D lhs, const double& rhs) { return lhs -= rhs; }
inline Vector2D operator-(const double& lhs, Vector2D rhs) { return rhs -= lhs; }
inline Vector2D operator-(Vector2D rhs) { return rhs *= -1; }

inline Vector2D operator*(Vector2D lhs, const double& rhs) { return lhs *= rhs; }
inline Vector2D operator*(const double& lhs, Vector2D rhs) { return rhs *= lhs; }

inline Vector2D operator/(Vector2D lhs, const double& rhs) { return lhs /= rhs; }

inline bool operator==(const Vector2D& lhs, const Vector2D& rhs) { return lhs.X == rhs.X && lhs.Y == rhs.Y; }
inline bool operator!=(const Vector2D& lhs, const Vector2D& rhs) { return !operator==(lhs, rhs); }