#include <iostream>
#include <string>

#include "ColorD.h"
#include "Vector3D.h"
#include "Matrix3x3D.h"
#include "Vertex.h"
#include "Triangle.h"
#include "Ray.h"
#include "Intersections.h"

int main()
{
	std::cout << "Test!" << std::endl;
	std::string s;

	/*Matrix3x3D mat1(3, -1, 2, 0, 4, 1, 5, 3, -4), mat2(6, 11, 4, 9, 0, 3, 1, 6, 2);
	mat.Invert();
	Matrix3x3D mat = mat1 * mat2 * Matrix3x3D::Identity();
	std::cout << mat.Determinant() << std::endl;
	std::cout << mat[0][0] << ", " << mat[0][1] << ", " << mat[0][2] << std::endl;
	std::cout << mat[1][0] << ", " << mat[1][1] << ", " << mat[1][2] << std::endl;
	std::cout << mat[2][0] << ", " << mat[2][1] << ", " << mat[2][2] << std::endl << std::endl;

	Vector3D vec = mat * Vector3D(4, 5, 6);*/

	/*Vector3D in(0.707107, -0.707107, 0), n(0, 1, 0);
	in.Refract(n, 0.9, 1);

	std::cout << in.X << ", " << in.Y << ", " << in.Z << std::endl;*/

	Triangle tri(Vertex(Vector3D(0, 0, 0), Vector3D(0, 0, 1), ColorD(1.0, 1.0, 1.0), Vector3D(0, 0, 0)),
				 Vertex(Vector3D(1, 0, 0), Vector3D(0, 0, 1), ColorD(1.0, 1.0, 1.0), Vector3D(0, 0, 0)),
				 Vertex(Vector3D(0, 2, 0), Vector3D(0, 0, 1), ColorD(1.0, 1.0, 1.0), Vector3D(0, 0, 0)));
	BoundingBox aabb(Vector3D(2.01, 0, 0), Vector3D(1, 1, 1));

	std::cout << Intersects(tri, aabb) << std::endl;

	getline(std::cin, s);
}