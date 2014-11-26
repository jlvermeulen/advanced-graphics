#include <math.h>
#include "Intersections.h"

#define EPSILON 0.000001

// Möller-Trumbore intersection algorithm
bool Intersects(const RayD& ray, const TriangleD& triangle, double& t)
{
	Vector3D e1 = triangle.Vertices[1].Position - triangle.Vertices[0].Position;
	Vector3D e2 = triangle.Vertices[2].Position - triangle.Vertices[0].Position;
	Vector3D p = ray.Direction.Cross(e2);
	double det = e1.Dot(p);

	if (det > - EPSILON && det < EPSILON)
		return false;

	double invDet = 1 / det;
	Vector3D r = ray.Origin - triangle.Vertices[0].Position;
	double u = r.Dot(p) * invDet;

	if (u < 0 || u > 1)
		return false;

	Vector3D q = r.Cross(e1);
	double v = ray.Direction.Dot(q) * invDet;

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

// based on "Fast 3D Triangle-Box Overlap Testing", http://fileadmin.cs.lth.se/cs/Personal/Tomas_Akenine-Moller/code/tribox3.txt
/*======================== X-tests ========================*/

#define AXISTEST_X01(a, b, fa, fb)										\
	p0 = a * v0.Y - b * v0.Z;											\
	p2 = a * v2.Y - b * v2.Z;											\
	if(p0 < p2) { min = p0; max = p2; } else { min = p2; max = p0; }	\
	rad = fa * boundingBox.Halfsize.Y + fb * boundingBox.Halfsize.Z;	\
	if(min > rad || max < - rad) return false;

#define AXISTEST_X2(a, b, fa, fb)										\
	p0 = a * v0.Y - b * v0.Z;											\
	p1 = a * v1.Y - b * v1.Z;											\
	if(p0 < p1) { min = p0; max = p1; } else { min = p1; max = p0; }	\
	rad = fa * boundingBox.Halfsize.Y + fb * boundingBox.Halfsize.Z;	\
	if(min > rad || max < - rad) return false;

/*======================== Y-tests ========================*/

#define AXISTEST_Y02(a, b, fa, fb)										\
	p0 = -a * v0.X + b * v0.Z;											\
	p2 = -a * v2.X + b * v2.Z;											\
	if(p0 < p2) { min = p0; max = p2; } else { min = p2; max = p0; }	\
	rad = fa * boundingBox.Halfsize.X + fb * boundingBox.Halfsize.Z;	\
	if(min > rad || max < - rad) return false;

#define AXISTEST_Y1(a, b, fa, fb)										\
	p0 = -a * v0.X + b * v0.Z;											\
	p1 = -a * v1.X + b * v1.Z;											\
	if(p0 < p1) { min = p0; max = p1; } else { min = p1; max = p0; }	\
	rad = fa * boundingBox.Halfsize.X + fb * boundingBox.Halfsize.Z;	\
	if(min > rad || max < - rad) return false;

/*======================== Z-tests ========================*/

#define AXISTEST_Z12(a, b, fa, fb)										\
	p1 = a * v1.X - b * v1.Y;											\
	p2 = a * v2.X - b * v2.Y;											\
	if(p2 < p1) { min = p2; max = p1; } else { min = p1; max = p2; }	\
	rad = fa * boundingBox.Halfsize.X + fb * boundingBox.Halfsize.Y;	\
	if(min > rad || max < - rad) return false;

#define AXISTEST_Z0(a, b, fa, fb)										\
	p0 = a * v0.X - b * v0.Y;											\
	p1 = a * v1.X - b * v1.Y;											\
	if(p0 < p1) { min = p0; max = p1; } else { min = p1; max = p0; }	\
	rad = fa * boundingBox.Halfsize.X + fb * boundingBox.Halfsize.Y;	\
	if(min > rad || max < - rad) return false;

#define FINDMINMAX(x0, x1, x2)			\
	min = max = x0;						\
	if(x1 < min) min = x1;				\
	if(x1 > max) max = x1;				\
	if(x2 < min) min = x2;				\
	if(x2 > max) max = x2;

bool Intersects(const TriangleD& triangle, const BoundingBox& boundingBox)
{
	Vector3D v0 = triangle.Vertices[0].Position - boundingBox.Center;
	Vector3D v1 = triangle.Vertices[1].Position - boundingBox.Center;
	Vector3D v2 = triangle.Vertices[2].Position - boundingBox.Center;

	Vector3D e0 = triangle.Vertices[1].Position - triangle.Vertices[0].Position;
	Vector3D e1 = triangle.Vertices[2].Position - triangle.Vertices[1].Position;
	Vector3D e2 = triangle.Vertices[0].Position - triangle.Vertices[2].Position;

	double p0, p1, p2, min, max, rad, fex, fey, fez;

	fex = fabsf(e0.X);
	fey = fabsf(e0.Y);
	fez = fabsf(e0.Z);
	AXISTEST_X01(e0.Z, e0.Y, fez, fey);
	AXISTEST_Y02(e0.Z, e0.X, fez, fex);
	AXISTEST_Z12(e0.Y, e0.X, fey, fex);

	fex = fabsf(e1.X);
	fey = fabsf(e1.Y);
	fez = fabsf(e1.Z);
	AXISTEST_X01(e1.Z, e1.Y, fez, fey);
	AXISTEST_Y02(e1.Z, e1.X, fez, fex);
	AXISTEST_Z0(e1.Y, e1.X, fey, fex);

	fex = fabsf(e2.X);
	fey = fabsf(e2.Y);
	fez = fabsf(e2.Z);
	AXISTEST_X2(e2.Z, e2.Y, fez, fey);
	AXISTEST_Y1(e2.Z, e2.X, fez, fex);
	AXISTEST_Z12(e2.Y, e2.X, fey, fex);

	FINDMINMAX(v0.X, v1.X, v2.X);
	if (min > boundingBox.Halfsize.X || max < -boundingBox.Halfsize.X)
		return false;

	FINDMINMAX(v0.Y, v1.Y, v2.Y);
	if (min > boundingBox.Halfsize.Y || max < -boundingBox.Halfsize.Y)
		return false;

	FINDMINMAX(v0.Z, v1.Z, v2.Z);
	if (min > boundingBox.Halfsize.Z || max < -boundingBox.Halfsize.Z)
		return false;

	Vector3D normal = e0.Cross(e1), vMin, vMax;
	if (normal.X > 0)
	{
		vMin.X = -boundingBox.Halfsize.X - v0.X;
		vMax.X = boundingBox.Halfsize.X - v0.X;
	}
	else
	{
		vMin.X = boundingBox.Halfsize.X - v0.X;
		vMax.X = -boundingBox.Halfsize.X - v0.X;
	}

	if (normal.Y > 0)
	{
		vMin.Y = -boundingBox.Halfsize.Y - v0.Y;
		vMax.Y = boundingBox.Halfsize.Y - v0.Y;
	}
	else
	{
		vMin.Y = boundingBox.Halfsize.Y - v0.Y;
		vMax.Y = -boundingBox.Halfsize.Y - v0.Y;
	}

	if (normal.Z > 0)
	{
		vMin.Z = -boundingBox.Halfsize.Z - v0.Z;
		vMax.Z = boundingBox.Halfsize.Z - v0.Z;
	}
	else
	{
		vMin.Z = boundingBox.Halfsize.Z - v0.Z;
		vMax.Z = -boundingBox.Halfsize.Z - v0.Z;
	}

	if (normal.Dot(vMin) > 0 || normal.Dot(vMax) < 0)
		return false;

	return true;
}