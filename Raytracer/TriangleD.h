#pragma once

#include "Vertex.h"

struct TriangleD
{
public:
	Vertex Vertices[3];

	TriangleD();
	TriangleD(Vertex v1, Vertex v2, Vertex v3);
	TriangleD(Vertex v[3]);

private:
};