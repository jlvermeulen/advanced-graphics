#pragma once

#include "Color3F.h"
#include "Vertex.h"

struct Triangle
{
private:
	float d00, d01, d11, invDenom;
	Vector3F v0, v1;

public:
	Vertex Vertices[3];
	float Area;

	Triangle() { }

	Triangle(const Vertex& v1, const Vertex& v2, const Vertex& v3)
	{
		Vertices[0] = v1;
		Vertices[1] = v2;
		Vertices[2] = v3;

		PreCalc();
	}

	Triangle(const Vertex v[3])
	{
		Vertices[0] = v[0];
		Vertices[1] = v[1];
		Vertices[2] = v[2];

		PreCalc();
	}

	void PreCalc()
	{
		v0 = Vertices[1].Position - Vertices[0].Position;
		v1 = Vertices[2].Position - Vertices[0].Position;

		Area = 0.5f * Vector3F::Cross(v0, v1).Length();
		
		d00 = Vector3F::Dot(v0, v0);
		d01 = Vector3F::Dot(v0, v1);
		d11 = Vector3F::Dot(v1, v1);
		invDenom = 1.0f / (d00 * d11 - d01 * d01);
	}

	Vector3F Interpolate(Vector3F point) const
	{
		Vector3F v2 = point - Vertices[0].Position;
		float d20 = Vector3F::Dot(v2, v0);
		float d21 = Vector3F::Dot(v2, v1);

		float a2 = (d11 * d20 - d01 * d21) * invDenom;
		float a3 = (d00 * d21 - d01 * d20) * invDenom;
		float a1 = 1.0f - a2 - a3;

		return Vector3F(a1, a2, a3);
	}
  
	Vector3F surfaceNormal(Vector3F point) const
	{
		Vector3F factors = Interpolate(point);
		return factors.X * Vertices[0].Normal + factors.Y * Vertices[1].Normal + factors.Z * Vertices[2].Normal;
	}
};

inline bool operator==(const Triangle& lhs, const Triangle& rhs) { return lhs.Vertices[0].Position == rhs.Vertices[0].Position &&
																		  lhs.Vertices[1].Position == rhs.Vertices[1].Position &&
																		  lhs.Vertices[2].Position == rhs.Vertices[2].Position; }

inline bool operator!=(const Triangle& lhs, const Triangle& rhs) { return !operator==(lhs, rhs); }
