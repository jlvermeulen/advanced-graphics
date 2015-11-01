#include "BoundingSphere.h"

#include <algorithm>
#include <math.h>

float dist(Vector3F v1, Vector3F v2){ return (v1 - v2).Length(); }

BoundingSphere BoundingSphere::FromTriangles(const std::vector<Triangle*>& triangles)
{
	// center
	Vector3F center;
	for (Triangle* t : triangles)
		center += (t->Vertices[0].Position + t->Vertices[1].Position + t->Vertices[2].Position);
	center /= triangles.size() * 3;

	//radius
	float radius = 0;
	for (Triangle* t : triangles)
	{
		float m1 = std::max(dist(t->Vertices[1].Position, center), dist(t->Vertices[2].Position, center));
		radius = std::max(radius, std::max(dist(t->Vertices[0].Position, center), m1));
	}

	return BoundingSphere(center, radius);
}