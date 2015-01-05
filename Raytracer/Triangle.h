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

	Vector3D Interpolate(Vector3D point) const
	{
		Vector3D p1 = Vertices[0].Position - point;
		Vector3D p2 = Vertices[1].Position - point;
		Vector3D p3 = Vertices[2].Position - point;

		double a = Vector3D::Cross(Vertices[0].Position - Vertices[1].Position, Vertices[0].Position - Vertices[2].Position).Length();
		double aInv = 1 / a;
		double a1 = Vector3D::Cross(p2, p3).Length() * aInv;
		double a2 = Vector3D::Cross(p3, p1).Length() * aInv;
		double a3 = Vector3D::Cross(p1, p2).Length() * aInv;

		return Vector3D(a1, a2, a3);
	}

	ColorD surfaceColor(Vector3D point) const
	{
		Vector3D factors = Interpolate(point);

		return factors.X * Vertices[0].Color + factors.Y * Vertices[1].Color + factors.Z * Vertices[2].Color;
	}
  
	Vector3D surfaceNormal(Vector3D point) const
	{
		Vector3D factors = Interpolate(point);

		return factors.X * Vertices[0].Normal + factors.Y * Vertices[1].Normal + factors.Z * Vertices[2].Normal;
	}

private:
};