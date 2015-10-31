#include "Vector2F.h"

Vector2F::Vector2F() : X(0), Y(0) { }
Vector2F::Vector2F(const float& x, const float& y) : X(x), Y(y) { }
Vector2F::Vector2F(const float vec[]) { X = vec[0]; Y = vec[1]; }

float Vector2F::Length() { return sqrt(LengthSquared()); }
float Vector2F::LengthSquared() { return Dot(*this, *this); }

float Vector2F::Dot(const Vector2F& other) const { return Dot(*this, other); }
float Vector2F::Dot(const Vector2F& lhs, const Vector2F& rhs) { return lhs.X * rhs.X + lhs.Y * rhs.Y; }

Vector2F& Vector2F::operator+=(const Vector2F& rhs)
{
	X += rhs.X;
	Y += rhs.Y;
	return *this;
}

Vector2F& Vector2F::operator+=(const float& rhs)
{
	X += rhs;
	Y += rhs;
	return *this;
}

Vector2F& Vector2F::operator-=(const Vector2F& rhs)
{
	X -= rhs.X;
	Y -= rhs.Y;
	return *this;
}

Vector2F& Vector2F::operator-=(const float& rhs)
{
	X -= rhs;
	Y -= rhs;
	return *this;
}

Vector2F& Vector2F::operator*=(const float& rhs)
{
	X *= rhs;
	Y *= rhs;
	return *this;
}

Vector2F& Vector2F::operator/=(const float& rhs)
{
	X /= rhs;
	Y /= rhs;
	return *this;
}

float& Vector2F::operator[](int index) { return index == 0 ? X : Y; }
const float& Vector2F::operator[](int index) const { return index == 0 ? X : Y; };
