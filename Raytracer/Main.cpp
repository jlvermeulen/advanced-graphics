#include <iostream>
#include <string>
#include "Vector3D.h"
#include "TriangleD.h"
#include "Matrix3x3D.h"

using namespace std;

int main()
{
	cout << "Test!" << endl;
	string s;
	/*int x, y, z;
	Vector3D v1, v2, v3;
	TriangleD t;
	getline(cin, s);
	x = atoi(s.c_str());
	getline(cin, s);
	y = atoi(s.c_str());
	getline(cin, s);
	z = atoi(s.c_str());
	v1 = Vector3D(x, y, z);
	getline(cin, s);
	x = atoi(s.c_str());
	getline(cin, s);
	y = atoi(s.c_str());
	getline(cin, s);
	z = atoi(s.c_str());
	v2 = Vector3D(x, y, z);
	getline(cin, s);
	x = atoi(s.c_str());
	getline(cin, s);
	y = atoi(s.c_str());
	getline(cin, s);
	z = atoi(s.c_str());
	v3 = Vector3D(x, y, z);

	v1 += v2;
	v2 -= v3;
	v3 = v3.Cross(v1, v2);

	Vector3D a[3] = { v1, v2, v3 };
	t = TriangleD(a);
	for (int i = 0; i < 3; i++)
		cout << "X: " << t.Vertices[i].X << ", Y: " << t.Vertices[i].Y << ", Z: " << t.Vertices[i].Z << endl;

	cout << v1.Dot(v2) << endl;*/

	Matrix3x3D mat(1, 2, 3, 4, 5, 6, 7, 8, 8);
	mat.Invert();
	cout << mat.Determinant() << endl;
	cout << mat[0][0] << ", " << mat[0][1] << ", " << mat[0][2] << endl;
	cout << mat[1][0] << ", " << mat[1][1] << ", " << mat[1][2] << endl;
	cout << mat[2][0] << ", " << mat[2][1] << ", " << mat[2][2] << endl << endl;

	Vector3D vec = mat * Vector3D(4, 5, 6);
	cout << vec.X << ", " << vec.Y << ", " << vec.Z << endl;

	getline(cin, s);
}