#include <algorithm>
#include "Color3F.h"

Color3F::Color3F() : R(0), G(0), B(0) { }
Color3F::Color3F(float r, float g, float b) : R(r), G(g), B(b) { }

Color3F::Color3F(unsigned char r, unsigned char g, unsigned char b)
{
  R = (float) r / 255;
  G = (float) g / 255;
  B = (float) b / 255;
}

void Color3F::Clip()
{
	R = std::max(0.0f, std::min(1.0f, R));
	G = std::max(0.0f, std::min(1.0f, G));
	B = std::max(0.0f, std::min(1.0f, B));
}

float Color3F::Magnitude()
{
	return std::max(R, std::max(G, B));
}

Color3F Color3F::Mix(const Color3F& c1, const Color3F& c2, const float& w1, const float& w2) { return w1 * c1 + w2 * c2; }

Color3F& Color3F::operator+=(const Color3F& rhs)
{
	R += rhs.R;
	G += rhs.G;
	B += rhs.B;
	return *this;
}

Color3F& Color3F::operator+=(const float& rhs)
{
	R += rhs;
	G += rhs;
	B += rhs;
	return *this;
}

Color3F& Color3F::operator*=(const Color3F& rhs)
{
	R *= rhs.R;
	G *= rhs.G;
	B *= rhs.B;
	return *this;
}

Color3F& Color3F::operator*=(const float& rhs)
{
	R *= rhs;
	G *= rhs;
	B *= rhs;
	return *this;
}

Color3F& Color3F::operator/=(const Color3F& rhs)
{
	R /= rhs.R;
	G /= rhs.G;
	B /= rhs.B;
	return *this;
}

Color3F& Color3F::operator/=(const float& rhs)
{
	R /= rhs;
	G /= rhs;
	B /= rhs;
	return *this;
}

float& Color3F::operator[](int index) { return index == 0 ? R : index == 1 ? G : B; }
const float& Color3F::operator[](int index) const { return index == 0 ? R : index == 1 ? G : B; };
