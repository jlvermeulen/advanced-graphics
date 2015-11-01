#pragma once

class Matrix3x3F
{
public:
	Matrix3x3F();
	Matrix3x3F(const float elems[3][3]);
	Matrix3x3F(const float elems[3]);
	Matrix3x3F(const float& m00, const float& m01, const float& m02,
			   const float& m10, const float& m11, const float& m12,
			   const float& m20, const float& m21, const float& m22);

	static Matrix3x3F Identity();
	static Matrix3x3F CreateRotationX(float angle);
	static Matrix3x3F CreateRotationY(float angle);
	static Matrix3x3F CreateRotationZ(float angle);

	void Transpose();
	Matrix3x3F Transposed() const;
	static Matrix3x3F Transpose(const Matrix3x3F& matrix);

	float Determinant();

	void Invert();
	Matrix3x3F Inverse() const;
	static Matrix3x3F Invert(const Matrix3x3F& matrix);

	float* operator[](int i);
	const float* operator[](int i) const;

	Matrix3x3F& operator*=(const Matrix3x3F& rhs);

private:
	float elems[3][3];
};

inline Matrix3x3F operator*(Matrix3x3F lhs, const Matrix3x3F& rhs) { return lhs *= rhs; }
