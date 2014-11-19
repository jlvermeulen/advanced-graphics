#pragma once

#include <algorithm>
#include "Color.h"

Color::Color() : R(0), G(0), B(0) { }
Color::Color(double r, double g, double b) : R(r), G(g), B(b) { }

Color::Color(unsigned char r, unsigned char g, unsigned char b)
{
	R = (double)r / 255;
	G = (double)g / 255;
	B = (double)b / 255;
}

void Color::Clip()
{
	R = std::max(0.0, std::min(1.0, R));
	G = std::max(0.0, std::min(1.0, G));
	B = std::max(0.0, std::min(1.0, B));
}

Color Color::Mix(Color c1, Color c2, const double& w1, const double& w2) { return w1 * c1 + w2 * c2; }

Color& Color::operator+=(const Color& rhs)
{
	R += rhs.R;
	G += rhs.G;
	B += rhs.B;
	return *this;
}

Color& Color::operator+=(const double& rhs)
{
	R += rhs;
	G += rhs;
	B += rhs;
	return *this;
}

Color& Color::operator*=(const Color& rhs)
{
	R *= rhs.R;
	G *= rhs.G;
	B *= rhs.B;
	return *this;
}

Color& Color::operator*=(const double& rhs)
{
	R *= rhs;
	G *= rhs;
	B *= rhs;
	return *this;
}

Color& Color::operator/=(const Color& rhs)
{
	R /= rhs.R;
	G /= rhs.G;
	B /= rhs.B;
	return *this;
}

Color& Color::operator/=(const double& rhs)
{
	R /= rhs;
	G /= rhs;
	B /= rhs;
	return *this;
}