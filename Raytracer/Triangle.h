#pragma once

#include "ColorD.h"
#include "Vertex.h"

struct Triangle
{
public:
	Vertex Vertices[3];

	Triangle() { }

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

  ColorD surfaceColor(Vector3D point) const
  {
    double d0 = (Vertices[0].Position - point).Length();
    double d1 = (Vertices[1].Position - point).Length();
    double d2 = (Vertices[2].Position - point).Length();

    return (d0 * Vertices[0].Color + d1 * Vertices[1].Color + d2 * Vertices[2].Color) / (d0 + d1 + d2);
  }
  
  Vector3D surfaceNormal(Vector3D point) const
  {
    double d0 = (Vertices[0].Position - point).Length();
    double d1 = (Vertices[1].Position - point).Length();
    double d2 = (Vertices[2].Position - point).Length();

    return (d0 * Vertices[0].Normal + d1 * Vertices[1].Normal + d2 * Vertices[2].Normal) / (d0 + d1 + d2);
  }

private:
};