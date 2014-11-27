#pragma once

class Matrix3x3D
{
public:
	Matrix3x3D();
	Matrix3x3D(const double elems[3][3]);
	Matrix3x3D(const double elems[3]);
	Matrix3x3D(const double& m00, const double& m01, const double& m02,
			   const double& m10, const double& m11, const double& m12,
			   const double& m20, const double& m21, const double& m22);

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