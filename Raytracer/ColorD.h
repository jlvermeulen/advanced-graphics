#pragma once

#define MIN_INTENSITY 0.001

class ColorD
{
public:
	double R, G, B;

	ColorD();
	ColorD(double r, double g, double b);
	ColorD(unsigned char r, unsigned char g, unsigned char b);

	inline bool ColorD::IsSignificant() const
	{
		return (R >= MIN_INTENSITY ||
				G >= MIN_INTENSITY ||
				B >= MIN_INTENSITY);
	}

	void Clip();

	static ColorD Mix(const ColorD& c1, const ColorD& c2, const double& w1, const double& w2);

	ColorD& operator+=(const ColorD& rhs);
	ColorD& operator+=(const double& rhs);

	ColorD& operator*=(const ColorD& rhs);
	ColorD& operator*=(const double& rhs);

	ColorD& operator/=(const ColorD& rhs);
	ColorD& operator/=(const double& rhs);

	double& operator[](int index);
	const double& operator[](int index) const;

private:
};

inline ColorD operator+(ColorD lhs, const ColorD& rhs) { return lhs += rhs; }
inline ColorD operator+(ColorD lhs, const double& rhs) { return lhs += rhs; }
inline ColorD operator+(const double& lhs, ColorD rhs) { return rhs += lhs; }

inline ColorD operator*(ColorD lhs, const ColorD& rhs) { return lhs *= rhs; }
inline ColorD operator*(ColorD lhs, const double& rhs) { return lhs *= rhs; }
inline ColorD operator*(const double& lhs, ColorD rhs) { return rhs *= lhs; }

inline ColorD operator/(ColorD lhs, const ColorD& rhs) { return lhs /= rhs; }
inline ColorD operator/(ColorD lhs, const double& rhs) { return lhs /= rhs; }