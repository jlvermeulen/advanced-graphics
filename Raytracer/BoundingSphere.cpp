#include "BoundingSphere.h"

#include <algorithm>
#include <math.h>

BoundingSphere::BoundingSphere()
{
}

using std::max;
BoundingSphere::~BoundingSphere()
{
}
BoundingSphere BoundingSphere::FromTriangles(const std::deque<Triangle>& triangles)
{
	// center
	Vector3D center = Vector3D(0.0, 0.0, 0.0);
	for (Triangle t : triangles)
		center += (t.Vertices[0].Position + t.Vertices[1].Position + t.Vertices[2].Position);
	center /= triangles.size()*3;
	//radius
	double radius = 0;
	for (Triangle t : triangles)
		radius = max(radius, max(dist(t.Vertices[0].Position, center), max(dist(t.Vertices[1].Position, center), dist(t.Vertices[2].Position, center))));
	return BoundingSphere(center, radius);
}

double dist(Vector3D v1, Vector3D v2){return (v1 - v2).Length;}