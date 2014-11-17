#pragma once

#include "Vector3D.h"

struct TriangleD
{
public:
	TriangleD() { }

	TriangleD(Vector3D v1, Vector3D v2, Vector3D v3)
	{
		Vertices[0] = v1;
		Vertices[1] = v2;
		Vertices[2] = v3;
	}

	TriangleD(Vector3D v[3])
	{
		Vertices[0] = v[0];
		Vertices[1] = v[1];
		Vertices[2] = v[2];
	}

	~TriangleD() { }
	Vector3D Vertices[3];

private:
};