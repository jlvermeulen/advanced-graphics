#pragma once

struct Color
{
public:
	double R, G, B;

	Color() { R = G = B = 0; }
	Color(double r, double g, double b) : R(r), G(g), B(b) { }

	Color(unsigned char r, unsigned char g, unsigned char b)
	{
		R = (double)r / 255;
		G = (double)g / 255;
		B = (double)b / 255;
	}

	static Color Mix(const Color& c1, const Color& c2, const double& w1, const double& w2) {}
private:
};

inline Color& operator*(const Color& lhs, const Color& rhs) {}
inline Color& operator+(const Color& lhs, const Color& rhs) {}