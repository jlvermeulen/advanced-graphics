#pragma once

#define MIN_INTENSITY 0.001f


class Color3F
{
public:
	float R, G, B;

	Color3F();
	Color3F(float r, float g, float b);
	Color3F(unsigned char r, unsigned char g, unsigned char b);

	inline bool Color3F::IsSignificant() const
	{
		return (R >= MIN_INTENSITY ||
				G >= MIN_INTENSITY ||
				B >= MIN_INTENSITY);
	}

	void Clip();
	float Magnitude();

	static Color3F Mix(const Color3F& c1, const Color3F& c2, const float& w1, const float& w2);

	Color3F& operator+=(const Color3F& rhs);
	Color3F& operator+=(const float& rhs);

	Color3F& operator*=(const Color3F& rhs);
	Color3F& operator*=(const float& rhs);

	Color3F& operator/=(const Color3F& rhs);
	Color3F& operator/=(const float& rhs);

	float& operator[](int index);
	const float& operator[](int index) const;

private:
};

inline Color3F operator+(Color3F lhs, const Color3F& rhs) { return lhs += rhs; }
inline Color3F operator+(Color3F lhs, const float& rhs) { return lhs += rhs; }
inline Color3F operator+(const float& lhs, Color3F rhs) { return rhs += lhs; }

inline Color3F operator*(Color3F lhs, const Color3F& rhs) { return lhs *= rhs; }
inline Color3F operator*(Color3F lhs, const float& rhs) { return lhs *= rhs; }
inline Color3F operator*(const float& lhs, Color3F rhs) { return rhs *= lhs; }

inline Color3F operator/(Color3F lhs, const Color3F& rhs) { return lhs /= rhs; }
inline Color3F operator/(Color3F lhs, const float& rhs) { return lhs /= rhs; }
