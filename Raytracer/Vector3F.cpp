#include "Vector3F.h"

Vector3F::Vector3F() : X(0), Y(0), Z(0) { }
Vector3F::Vector3F(const float& x, const float& y, const float& z) : X(x), Y(y), Z(z) { }
Vector3F::Vector3F(const float vec[]) { X = vec[0]; Y = vec[1]; Z = vec[2]; }

float Vector3F::Length() const { return sqrtf(LengthSquared()); }
float Vector3F::LengthSquared() const { return Dot(*this, *this); }

float Vector3F::Dot(const Vector3F& other) const { return Dot(*this, other); }
float Vector3F::Dot(const Vector3F& lhs, const Vector3F& rhs) { return lhs.X * rhs.X + lhs.Y * rhs.Y + lhs.Z * rhs.Z; }

Vector3F Vector3F::Cross(const Vector3F& other) const { return Cross(*this, other); }
Vector3F Vector3F::Cross(const Vector3F& lhs, const Vector3F& rhs) { return Vector3F(lhs.Y * rhs.Z - lhs.Z * rhs.Y, lhs.Z * rhs.X - lhs.X * rhs.Z, lhs.X * rhs.Y - lhs.Y * rhs.X); }

void Vector3F::Normalise() { (*this) /= Length(); }
Vector3F Vector3F::Normalise(Vector3F vector) { vector.Normalise(); return vector; }

Vector3F& Vector3F::operator+=(const Vector3F& rhs)
{
	X += rhs.X;
	Y += rhs.Y;
	Z += rhs.Z;
	return *this;
}

Vector3F& Vector3F::operator+=(const float& rhs)
{
	X += rhs;
	Y += rhs;
	Z += rhs;
	return *this;
}

Vector3F& Vector3F::operator-=(const Vector3F& rhs)
{
	X -= rhs.X;
	Y -= rhs.Y;
	Z -= rhs.Z;
	return *this;
}

Vector3F& Vector3F::operator-=(const float& rhs)
{
	X -= rhs;
	Y -= rhs;
	Z -= rhs;
	return *this;
}

Vector3F& Vector3F::operator*=(const float& rhs)
{
	X *= rhs;
	Y *= rhs;
	Z *= rhs;
	return *this;
}

Vector3F& Vector3F::operator*=(const Matrix3x3F& lhs)
{
	Vector3F v(*this);
	X = v.Dot(Vector3F(lhs[0]));
	Y = v.Dot(Vector3F(lhs[1]));
	Z = v.Dot(Vector3F(lhs[2]));
	return *this;
}

Vector3F& Vector3F::operator/=(const float& rhs)
{
	X /= rhs;
	Y /= rhs;
	Z /= rhs;
	return *this;
}

float& Vector3F::operator[](int index) { return index == 0 ? X : index == 1 ? Y : Z; }
const float& Vector3F::operator[](int index) const { return index == 0 ? X : index == 1 ? Y : Z; };

Vector3F Vector3F::Forward = Vector3F(0, 0, -1);
Vector3F Vector3F::Right = Vector3F(1, 0, 0);
Vector3F Vector3F::Up = Vector3F(0, 1, 0);
