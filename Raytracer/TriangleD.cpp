#pragma once

#include "TriangleD.h"

TriangleD::TriangleD() { }

TriangleD::TriangleD(Vertex v1, Vertex v2, Vertex v3)
{
	Vertices[0] = v1;
	Vertices[1] = v2;
	Vertices[2] = v3;
}

TriangleD::TriangleD(Vertex v[3])
{
	Vertices[0] = v[0];
	Vertices[1] = v[1];
	Vertices[2] = v[2];
}