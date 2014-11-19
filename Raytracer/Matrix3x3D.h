#pragma once

struct Matrix3x3D
{
public:
	Matrix3x3D();
	Matrix3x3D(double elems[3][3]);
	Matrix3x3D(double elems[3]);
	Matrix3x3D(double m00, double m01, double m02, double m10, double m11, double m12, double m20, double m21, double m22);

	static Matrix3x3D Identity();

	void Transpose();
	Matrix3x3D Transposed() const;
	static Matrix3x3D Transpose(const Matrix3x3D& matrix);

	double Determinant();

	void Invert();
	Matrix3x3D Inverse() const;
	static Matrix3x3D Invert(const Matrix3x3D& matrix);

	double* operator[](int i);
	const double* operator[](int i) const;

	Matrix3x3D& operator*=(const Matrix3x3D& rhs);

private:
	double elems[3][3];
};

inline Matrix3x3D operator*(Matrix3x3D lhs, const Matrix3x3D& rhs) { return lhs *= rhs; }