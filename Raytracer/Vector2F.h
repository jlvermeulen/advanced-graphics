#pragma once

#include <math.h>

class Vector2F
{
public:
	float X, Y;

	Vector2F();
	Vector2F(const float& x, const float& y);
	Vector2F(const float vec[]);

	float Length();
  float LengthSquared();

	float Dot(const Vector2F& other) const;
	static float Dot(const Vector2F& lhs, const Vector2F& rhs);

	void Normalise();
	static Vector2F Normalise(Vector2F vector);

	Vector2F& operator+=(const Vector2F& rhs);
	Vector2F& operator+=(const float& rhs);

	Vector2F& operator-=(const Vector2F& rhs);
	Vector2F& operator-=(const float& rhs);
	
	Vector2F& operator*=(const float& rhs);

	Vector2F& operator/=(const float& rhs);

	float& operator[](int index);
	const float& operator[](int index) const;
};

inline Vector2F operator+(Vector2F lhs, const Vector2F& rhs) { return lhs += rhs; }
inline Vector2F operator+(Vector2F lhs, const float& rhs) { return lhs += rhs; }
inline Vector2F operator+(const float& lhs, Vector2F rhs) { return rhs += lhs; }

inline Vector2F operator-(Vector2F lhs, const Vector2F& rhs) { return lhs -= rhs; }
inline Vector2F operator-(Vector2F lhs, const float& rhs) { return lhs -= rhs; }
inline Vector2F operator-(const float& lhs, Vector2F rhs) { return rhs -= lhs; }
inline Vector2F operator-(Vector2F rhs) { return rhs *= -1; }

inline Vector2F operator*(Vector2F lhs, const float& rhs) { return lhs *= rhs; }
inline Vector2F operator*(const float& lhs, Vector2F rhs) { return rhs *= lhs; }

inline Vector2F operator/(Vector2F lhs, const float& rhs) { return lhs /= rhs; }

inline bool operator==(const Vector2F& lhs, const Vector2F& rhs) { return lhs.X == rhs.X && lhs.Y == rhs.Y; }
inline bool operator!=(const Vector2F& lhs, const Vector2F& rhs) { return !operator==(lhs, rhs); }
