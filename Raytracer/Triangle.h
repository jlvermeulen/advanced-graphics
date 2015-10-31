#pragma once

#include "Color3F.h"
#include "Vertex.h"

struct Triangle
{
public:
	Vertex Vertices[3];
	float Area;
	Vector3F Center;

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

	Vector3F Interpolate(Vector3F point) const
	{
		Vector3F p1 = Vertices[0].Position - point;
		Vector3F p2 = Vertices[1].Position - point;
		Vector3F p3 = Vertices[2].Position - point;

		float a = Vector3F::Cross(Vertices[0].Position - Vertices[1].Position, Vertices[0].Position - Vertices[2].Position).Length();
		float aInv = 1 / a;
		float a1 = Vector3F::Cross(p2, p3).Length() * aInv;
		float a2 = Vector3F::Cross(p3, p1).Length() * aInv;
		float a3 = Vector3F::Cross(p1, p2).Length() * aInv;

		return Vector3F(a1, a2, a3);
	}
  
	Vector3F surfaceNormal(Vector3F point) const
	{
		Vector3F factors = Interpolate(point);

		return factors.X * Vertices[0].Normal + factors.Y * Vertices[1].Normal + factors.Z * Vertices[2].Normal;
	}

	void CalculateArea() { Area = 0.5f * Vector3F::Cross(Vertices[1].Position - Vertices[0].Position, Vertices[2].Position - Vertices[0].Position).Length(); }
	void CalculateCenter() { Center = (Vertices[0].Position + Vertices[1].Position + Vertices[2].Position) / 3; }

private:
};

inline bool operator==(const Triangle& lhs, const Triangle& rhs) { return lhs.Vertices[0].Position == rhs.Vertices[0].Position &&
																		  lhs.Vertices[1].Position == rhs.Vertices[1].Position &&
																		  lhs.Vertices[2].Position == rhs.Vertices[2].Position; }

inline bool operator!=(const Triangle& lhs, const Triangle& rhs) { return !operator==(lhs, rhs); }
