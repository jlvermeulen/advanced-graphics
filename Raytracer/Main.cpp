#include <iostream>
#include <string>

#include "Color.h"
#include "Vector3D.h"
#include "Matrix3x3D.h"
#include "Vertex.h"
#include "TriangleD.h"

int main()
{
	std::cout << "Test!" << std::endl;
	std::string s;
	/*int x, y, z;
	Vector3D v1, v2, v3;
	TriangleD t;
	getline(std::cin, s);
	x = atoi(s.c_str());
	getline(std::cin, s);
	y = atoi(s.c_str());
	getline(std::cin, s);
	z = atoi(s.c_str());
	v1 = Vector3D(x, y, z);
	getline(std::cin, s);
	x = atoi(s.c_str());
	getline(std::cin, s);
	y = atoi(s.c_str());
	getline(std::cin, s);
	z = atoi(s.c_str());
	v2 = Vector3D(x, y, z);
	getline(std::cin, s);
	x = atoi(s.c_str());
	getline(std::cin, s);
	y = atoi(s.c_str());
	getline(std::cin, s);
	z = atoi(s.c_str());
	v3 = Vector3D(x, y, z);

	v1 += v2;
	v2 -= v3;
	v3 = v3.Cross(v1, v2);

	Vector3D a[3] = { v1, v2, v3 };
	t = TriangleD(a);
	for (int i = 0; i < 3; i++)
		std::cout << "X: " << t.Vertices[i].X << ", Y: " << t.Vertices[i].Y << ", Z: " << t.Vertices[i].Z << std::endl;

	std::cout << v1.Dot(v2) << std::endl;*/

	Matrix3x3D mat1(3, -1, 2, 0, 4, 1, 5, 3, -4), mat2(6, 11, 4, 9, 0, 3, 1, 6, 2);
	//mat.Invert();
	Matrix3x3D mat = mat1 * mat2 * Matrix3x3D::Identity();
	std::cout << mat.Determinant() << std::endl;
	std::cout << mat[0][0] << ", " << mat[0][1] << ", " << mat[0][2] << std::endl;
	std::cout << mat[1][0] << ", " << mat[1][1] << ", " << mat[1][2] << std::endl;
	std::cout << mat[2][0] << ", " << mat[2][1] << ", " << mat[2][2] << std::endl << std::endl;

	Vector3D vec = mat * Vector3D(4, 5, 6);
	std::cout << vec.X << ", " << vec.Y << ", " << vec.Z << std::endl;

	getline(std::cin, s);
}