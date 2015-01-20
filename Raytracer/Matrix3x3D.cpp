#include "Matrix3x3D.h"
#include "Vector3D.h"

Matrix3x3D::Matrix3x3D()
{
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			this->elems[i][j] = 0;
}

Matrix3x3D::Matrix3x3D(const double elems[3][3])
{
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			this->elems[i][j] = elems[i][j];
}

Matrix3x3D::Matrix3x3D(const double elems[3])
{
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			this->elems[i][j] = elems[i * 3 + j];
}

Matrix3x3D::Matrix3x3D(const double& m00, const double& m01, const double& m02, const double& m10, const double& m11, const double& m12, const double& m20, const double& m21, const double& m22)
{
	elems[0][0] = m00;
	elems[0][1] = m01;
	elems[0][2] = m02;
	elems[1][0] = m10;
	elems[1][1] = m11;
	elems[1][2] = m12;
	elems[2][0] = m20;
	elems[2][1] = m21;
	elems[2][2] = m22;
}

Matrix3x3D Matrix3x3D::Identity() { return Matrix3x3D(1, 0, 0, 0, 1, 0, 0, 0, 1); }

void Matrix3x3D::Transpose()
{
	Matrix3x3D mat(*this);
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			elems[i][j] = mat[j][i];
}

Matrix3x3D Matrix3x3D::Transposed() const
{
	Matrix3x3D mat(*this);
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			mat[i][j] = elems[j][i];
	return mat;
}

Matrix3x3D Matrix3x3D::Transpose(const Matrix3x3D& matrix) { return matrix.Transposed(); }

double Matrix3x3D::Determinant()
{
	return elems[0][0] * elems[1][1] * elems[2][2] +
		   elems[1][0] * elems[2][1] * elems[0][2] +
		   elems[2][0] * elems[0][1] * elems[1][2] -
		   elems[0][0] * elems[2][1] * elems[1][2] -
		   elems[2][0] * elems[1][1] * elems[0][2] -
		   elems[1][0] * elems[0][1] * elems[2][2];
}

void Matrix3x3D::Invert()
{
	double det = this->Determinant();
	if (!det)
		return;

	Matrix3x3D copy(*this);

	elems[0][0] = (copy[1][1] * copy[2][2] - copy[1][2] * copy[2][1]) / det;
	elems[0][1] = (copy[0][2] * copy[2][1] - copy[0][1] * copy[2][2]) / det;
	elems[0][2] = (copy[0][1] * copy[1][2] - copy[0][2] * copy[1][1]) / det;
	elems[1][0] = (copy[1][2] * copy[2][0] - copy[1][0] * copy[2][2]) / det;
	elems[1][1] = (copy[0][0] * copy[2][2] - copy[0][2] * copy[2][0]) / det;
	elems[1][2] = (copy[0][2] * copy[1][0] - copy[0][0] * copy[1][2]) / det;
	elems[2][0] = (copy[1][0] * copy[2][1] - copy[1][1] * copy[2][0]) / det;
	elems[2][1] = (copy[0][1] * copy[2][0] - copy[0][0] * copy[2][1]) / det;
	elems[2][2] = (copy[0][0] * copy[1][1] - copy[0][1] * copy[1][0]) / det;
}

Matrix3x3D Matrix3x3D::Inverse() const { Matrix3x3D m(*this); m.Invert(); return m; }
Matrix3x3D Matrix3x3D::Invert(const Matrix3x3D& matrix) { return matrix.Inverse(); }

double* Matrix3x3D::operator[](int i) { return elems[i]; }
const double* Matrix3x3D::operator[](int i) const { return elems[i]; }

Matrix3x3D& Matrix3x3D::operator*=(const Matrix3x3D& rhs)
{
	Matrix3x3D mat(*this);
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			elems[i][j] = Vector3D::Dot(Vector3D(mat[i]), Vector3D(rhs[0][j], rhs[1][j], rhs[2][j]));

	return *this;
}

Matrix3x3D Matrix3x3D::CreateRotationX(double angle) { return Matrix3x3D(1, 0, 0, 0, cos(angle), -sin(angle), 0, sin(angle), cos(angle)); }
Matrix3x3D Matrix3x3D::CreateRotationY(double angle) { return Matrix3x3D(cos(angle), 0, sin(angle), 0, 1, 0, -sin(angle), 0, cos(angle)); }
Matrix3x3D Matrix3x3D::CreateRotationZ(double angle) { return Matrix3x3D(cos(angle), -sin(angle), 0, sin(angle), cos(angle), 0, 0, 0, 1); }