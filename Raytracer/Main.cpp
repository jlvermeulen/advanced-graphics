#include <iostream>
#include <string>
#include "Vector3D.h"
#include "Triangle.h"

using namespace std;

int main()
{
	cout << "Test!" << endl;
	string s;
	int x, y, z;
	Vector3D v1, v2, v3;
	Triangle t;
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
	t = Triangle(a);
	for (int i = 0; i < 3; i++)
		cout << "X: " << t.Vertices[i].X << ", Y: " << t.Vertices[i].Y << ", Z: " << t.Vertices[i].Z << endl;

	cout << v1.Dot(v2) << endl;

	getline(cin, s);
}