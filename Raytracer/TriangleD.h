#pragma once

#include "Vertex.h"

struct TriangleD
{
public:
	TriangleD() { }

	TriangleD(Vertex v1, Vertex v2, Vertex v3)
	{
		Vertices[0] = v1;
		Vertices[1] = v2;
		Vertices[2] = v3;
	}

	TriangleD(Vertex v[3])
	{
		Vertices[0] = v[0];
		Vertices[1] = v[1];
		Vertices[2] = v[2];
	}

	~TriangleD() { }
	Vertex Vertices[3];

private:
};