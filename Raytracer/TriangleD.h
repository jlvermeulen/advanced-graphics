#pragma once

#include "VertexD.h"

struct TriangleD
{
public:
	VertexD Vertices[3];

	TriangleD(const VertexD& v1, const VertexD& v2, const VertexD& v3)
	{
		Vertices[0] = v1;
		Vertices[1] = v2;
		Vertices[2] = v3;
	}

	TriangleD(const VertexD v[3])
	{
		Vertices[0] = v[0];
		Vertices[1] = v[1];
		Vertices[2] = v[2];
	}

private:
};