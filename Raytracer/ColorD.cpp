#include <algorithm>
#include "ColorD.h"

ColorD::ColorD() : R(0), G(0), B(0) { }
ColorD::ColorD(double r, double g, double b) : R(r), G(g), B(b) { }

ColorD::ColorD(unsigned char r, unsigned char g, unsigned char b)
{
	R = (double)r / 255;
	G = (double)g / 255;
	B = (double)b / 255;
}

void ColorD::Clip()
{
	R = std::max(0.0, std::min(1.0, R));
	G = std::max(0.0, std::min(1.0, G));
	B = std::max(0.0, std::min(1.0, B));
}

ColorD ColorD::Mix(const ColorD& c1, const ColorD& c2, const double& w1, const double& w2) { return w1 * c1 + w2 * c2; }

ColorD& ColorD::operator+=(const ColorD& rhs)
{
	R += rhs.R;
	G += rhs.G;
	B += rhs.B;
	return *this;
}

ColorD& ColorD::operator+=(const double& rhs)
{
	R += rhs;
	G += rhs;
	B += rhs;
	return *this;
}

ColorD& ColorD::operator*=(const ColorD& rhs)
{
	R *= rhs.R;
	G *= rhs.G;
	B *= rhs.B;
	return *this;
}

ColorD& ColorD::operator*=(const double& rhs)
{
	R *= rhs;
	G *= rhs;
	B *= rhs;
	return *this;
}

ColorD& ColorD::operator/=(const ColorD& rhs)
{
	R /= rhs.R;
	G /= rhs.G;
	B /= rhs.B;
	return *this;
}

ColorD& ColorD::operator/=(const double& rhs)
{
	R /= rhs;
	G /= rhs;
	B /= rhs;
	return *this;
}