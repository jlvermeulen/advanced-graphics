#include "Vector2D.h"

Vector2D::Vector2D() : X(0), Y(0) { }
Vector2D::Vector2D(const double& x, const double& y) : X(x), Y(y) { }
Vector2D::Vector2D(const double vec[]) { X = vec[0]; Y = vec[1]; }

double Vector2D::Length() { return sqrt(LengthSquared()); }
double Vector2D::LengthSquared() { return Dot(*this, *this); }

double Vector2D::Dot(const Vector2D& other) const { return Dot(*this, other); }
double Vector2D::Dot(const Vector2D& lhs, const Vector2D& rhs) { return lhs.X * rhs.X + lhs.Y * rhs.Y; }

Vector2D& Vector2D::operator+=(const Vector2D& rhs)
{
	X += rhs.X;
	Y += rhs.Y;
	return *this;
}

Vector2D& Vector2D::operator+=(const double& rhs)
{
	X += rhs;
	Y += rhs;
	return *this;
}

Vector2D& Vector2D::operator-=(const Vector2D& rhs)
{
	X -= rhs.X;
	Y -= rhs.Y;
	return *this;
}

Vector2D& Vector2D::operator-=(const double& rhs)
{
	X -= rhs;
	Y -= rhs;
	return *this;
}

Vector2D& Vector2D::operator*=(const double& rhs)
{
	X *= rhs;
	Y *= rhs;
	return *this;
}

Vector2D& Vector2D::operator/=(const double& rhs)
{
	X /= rhs;
	Y /= rhs;
	return *this;
}

double& Vector2D::operator[](int index) { return index == 0 ? X : Y; }
const double& Vector2D::operator[](int index) const { return index == 0 ? X : Y; };