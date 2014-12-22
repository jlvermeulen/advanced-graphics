#pragma once

#include "Vertex.h"

struct Triangle
{
public:
	Vertex Vertices[3];

	Triangle(const Vertex& v1, const Vertex& v2, const Vertex& v3)
	{
		Vertices[0] = v1;
		Vertices[1] = v2;
		Vertices[2] = v3;
	}

	Triangle(const Vertex v[3])
	{
		Vertices[0] = v[0];
		Vertices[1] = v[1];
		Vertices[2] = v[2];
	}

private:
};