#include "Matrix3x3F.h"
#include "Vector3F.h"

Matrix3x3F::Matrix3x3F()
{
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			this->elems[i][j] = 0;
}

Matrix3x3F::Matrix3x3F(const float elems[3][3])
{
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			this->elems[i][j] = elems[i][j];
}

Matrix3x3F::Matrix3x3F(const float elems[3])
{
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			this->elems[i][j] = elems[i * 3 + j];
}

Matrix3x3F::Matrix3x3F(const float& m00, const float& m01, const float& m02, const float& m10, const float& m11, const float& m12, const float& m20, const float& m21, const float& m22)
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

Matrix3x3F Matrix3x3F::Identity() { return Matrix3x3F(1, 0, 0, 0, 1, 0, 0, 0, 1); }

void Matrix3x3F::Transpose()
{
	Matrix3x3F mat(*this);
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			elems[i][j] = mat[j][i];
}

Matrix3x3F Matrix3x3F::Transposed() const
{
	Matrix3x3F mat(*this);
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			mat[i][j] = elems[j][i];
	return mat;
}

Matrix3x3F Matrix3x3F::Transpose(const Matrix3x3F& matrix) { return matrix.Transposed(); }

float Matrix3x3F::Determinant()
{
	return elems[0][0] * elems[1][1] * elems[2][2] +
		   elems[1][0] * elems[2][1] * elems[0][2] +
		   elems[2][0] * elems[0][1] * elems[1][2] -
		   elems[0][0] * elems[2][1] * elems[1][2] -
		   elems[2][0] * elems[1][1] * elems[0][2] -
		   elems[1][0] * elems[0][1] * elems[2][2];
}

void Matrix3x3F::Invert()
{
	float det = this->Determinant();
	if (!det)
		return;

	Matrix3x3F copy(*this);

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

Matrix3x3F Matrix3x3F::Inverse() const { Matrix3x3F m(*this); m.Invert(); return m; }
Matrix3x3F Matrix3x3F::Invert(const Matrix3x3F& matrix) { return matrix.Inverse(); }

float* Matrix3x3F::operator[](int i) { return elems[i]; }
const float* Matrix3x3F::operator[](int i) const { return elems[i]; }

Matrix3x3F& Matrix3x3F::operator*=(const Matrix3x3F& rhs)
{
	Matrix3x3F mat(*this);
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			elems[i][j] = Vector3F::Dot(Vector3F(mat[i]), Vector3F(rhs[0][j], rhs[1][j], rhs[2][j]));

	return *this;
}

Matrix3x3F Matrix3x3F::CreateRotationX(float angle) { return Matrix3x3F(1, 0, 0, 0, cosf(angle), -sinf(angle), 0, sinf(angle), cosf(angle)); }
Matrix3x3F Matrix3x3F::CreateRotationY(float angle) { return Matrix3x3F(cosf(angle), 0, sinf(angle), 0, 1, 0, -sinf(angle), 0, cosf(angle)); }
Matrix3x3F Matrix3x3F::CreateRotationZ(float angle) { return Matrix3x3F(cosf(angle), -sinf(angle), 0, sinf(angle), cosf(angle), 0, 0, 0, 1); }
