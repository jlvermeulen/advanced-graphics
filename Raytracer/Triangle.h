#pragma once

#include "ColorD.h"
#include "Vertex.h"
#include "BBox.h"
#include "Vector3.h"


struct Triangle
{
public:
	Vertex Vertices[3];
	double Area;
	Vector3D Center;
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

	BBox getBBox() const
	{
		Vector3 v0 = Vector3(Vertices[0].Position.X, Vertices[0].Position.Y, Vertices[0].Position.Z);
		Vector3 v1 = Vector3(Vertices[1].Position.X, Vertices[1].Position.Y, Vertices[1].Position.Z);
		Vector3 v2 = Vector3(Vertices[2].Position.X, Vertices[2].Position.Y, Vertices[2].Position.Z);
		const __declspec(align(16))Vector3 minv = min(min(v0, v1), v2);
		const __declspec(align(16))Vector3 maxv = max(max(v0, v1), v2);
		return BBox(minv, maxv);
	}
	Vector3 getCentroid()
	{
		CalculateCenter();
		return Vector3(Center.X, Center.Y, Center.Z);
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
  
	Vector3D surfaceNormal(Vector3D point) const
	{
		Vector3D factors = Interpolate(point);

		return factors.X * Vertices[0].Normal + factors.Y * Vertices[1].Normal + factors.Z * Vertices[2].Normal;
	}

	void CalculateArea() { Area = 0.5 * Vector3D::Cross(Vertices[1].Position - Vertices[0].Position, Vertices[2].Position - Vertices[0].Position).Length(); }
	void CalculateCenter() { Center = (Vertices[0].Position + Vertices[1].Position + Vertices[2].Position) / 3; }

private:
};

inline bool operator==(const Triangle& lhs, const Triangle& rhs) { return lhs.Vertices[0].Position == rhs.Vertices[0].Position &&
																		  lhs.Vertices[1].Position == rhs.Vertices[1].Position &&
																		  lhs.Vertices[2].Position == rhs.Vertices[2].Position; }

inline bool operator!=(const Triangle& lhs, const Triangle& rhs) { return !operator==(lhs, rhs); }
