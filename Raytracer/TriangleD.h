#pragma once

#include "Vertex.h"

struct TriangleD
{
public:
	Vertex Vertices[3];

	TriangleD(const Vertex& v1, const Vertex& v2, const Vertex& v3);
	TriangleD(const Vertex v[3]);

private:
};