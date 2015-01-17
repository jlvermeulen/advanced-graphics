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

	float area()
	{
		float x1, x2, x3, y1, y2, y3, z1, z2, z3;
		x1 = Vertices[0].Position.X;
		x2 = Vertices[1].Position.X;
		x3 = Vertices[2].Position.X;
		y1 = Vertices[0].Position.Y;
		y2 = Vertices[1].Position.Y;
		y3 = Vertices[2].Position.Y;
		z1 = Vertices[0].Position.Z;
		z2 = Vertices[1].Position.Z;
		z3 = Vertices[2].Position.Z;
		float t1 = (x2*y3 - x3*y2);
		float t2 = (x3*y1 - x1*y3);
		float t3 = (x1*y2 - x2*y1);
		return 0.5 * std::sqrtf(t1*t1 + t2*t2 + t3*t3);
	}

	Vector3D center() { return (Vertices[0].Position + Vertices[1].Position + Vertices[2].Position) / 3; }
private:
};