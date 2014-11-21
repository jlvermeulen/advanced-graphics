#pragma once

struct Color
{
public:
	double R, G, B;

	Color();
	Color(double r, double g, double b);
	Color(unsigned char r, unsigned char g, unsigned char b);

	void Clip();

	static Color Mix(const Color& c1, const Color& c2, const double& w1, const double& w2);

	Color& operator+=(const Color& rhs);
	Color& operator+=(const double& rhs);

	Color& operator*=(const Color& rhs);
	Color& operator*=(const double& rhs);

	Color& operator/=(const Color& rhs);
	Color& operator/=(const double& rhs);

private:
};

inline Color operator+(Color lhs, const Color& rhs) { return lhs += rhs; }
inline Color operator+(Color lhs, const double& rhs) { return lhs += rhs; }
inline Color operator+(const double& lhs, Color rhs) { return rhs += lhs; }

inline Color operator*(Color lhs, const Color& rhs) { return lhs *= rhs; }
inline Color operator*(Color lhs, const double& rhs) { return lhs *= rhs; }
inline Color operator*(const double& lhs, Color rhs) { return rhs *= lhs; }

inline Color operator/(Color lhs, const Color& rhs) { return lhs /= rhs; }
inline Color operator/(Color lhs, const double& rhs) { return lhs /= rhs; }