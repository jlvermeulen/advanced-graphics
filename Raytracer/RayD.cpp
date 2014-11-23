#include "RayD.h"

#define EPSILON 0.000001

RayD::RayD(const Vector3D& origin, const Vector3D& direction, const ColorD& color) : Origin(origin), Direction(direction), Color(color) { }

// Möller-Trumbore intersection algorithm
bool RayD::Intersects(const TriangleD& triangle, double& t) const
{
	Vector3D e1 = triangle.Vertices[1].Position - triangle.Vertices[0].Position;
	Vector3D e2 = triangle.Vertices[2].Position - triangle.Vertices[0].Position;
	Vector3D p = Direction.Cross(e2);
	double det = e1.Dot(p);

	if (det > -EPSILON && det < EPSILON)
		return false;

	double invDet = 1 / det;
	Vector3D r = Origin - triangle.Vertices[0].Position;
	double u = r.Dot(p) * invDet;

	if (u < 0 || u > 1)
		return false;

	Vector3D q = r.Cross(e1);
	double v = Direction.Dot(q) * invDet;

	if (v < 0 || u + v > 1)
		return false;

	double tt = e2.Dot(q) * invDet;

	if (tt > EPSILON)
	{
		t = tt;
		return true;
	}

	return false;
}