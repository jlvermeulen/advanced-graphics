//#include <math.h>
#include <algorithm>
#include <cmath>
#include "Intersections.h"

bool Intersects(const Ray& ray, const BoundingBox& bb, float& t)
{
	Vector3F lb = bb.Center - bb.Halfsize - ray.Origin, rt = bb.Center + bb.Halfsize - ray.Origin;

	float t1 = lb.X * ray.InverseDirection.X;
	float t2 = rt.X * ray.InverseDirection.X;

	float tmin = std::min(t1, t2);
	float tmax = std::max(t1, t2);

	t1 = lb.Y * ray.InverseDirection.Y;
	t2 = rt.Y * ray.InverseDirection.Y;

	tmin = std::max(tmin, std::min(t1, t2));
	tmax = std::min(tmax, std::max(t1, t2));

	t1 = lb.Z * ray.InverseDirection.Z;
	t2 = rt.Z * ray.InverseDirection.Z;

	tmin = std::max(tmin, std::min(t1, t2));
	tmax = std::min(tmax, std::max(t1, t2));

	if (tmax >= 0 && tmax >= tmin)
	{
		t = tmin;
		return true;
	}

	return false;
}

// M�ller-Trumbore intersection algorithm
bool Intersects(const Ray& ray, const Triangle& triangle, float& t)
{
	Vector3F e1 = triangle.Vertices[1].Position - triangle.Vertices[0].Position;
	Vector3F e2 = triangle.Vertices[2].Position - triangle.Vertices[0].Position;
	Vector3F p = ray.Direction.Cross(e2);
	float det = e1.Dot(p);

	if (det > - EPSILON && det < EPSILON)
		return false;

	float invDet = 1 / det;
	Vector3F r = ray.Origin - triangle.Vertices[0].Position;
	float u = r.Dot(p) * invDet;

	if (u < 0 || u > 1)
		return false;

	Vector3F q = r.Cross(e1);
	float v = ray.Direction.Dot(q) * invDet;

	if (v < 0 || u + v > 1)
		return false;

	float tt = e2.Dot(q) * invDet;

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

#define FINDMINMAX(x0, x1, x2)											\
	min = max = x0;														\
	if(x1 < min) min = x1;												\
	if(x1 > max) max = x1;												\
	if(x2 < min) min = x2;												\
	if(x2 > max) max = x2;

bool Intersects(const Triangle& triangle, const BoundingBox& boundingBox)
{
	Vector3F v0 = triangle.Vertices[0].Position - boundingBox.Center;
	Vector3F v1 = triangle.Vertices[1].Position - boundingBox.Center;
	Vector3F v2 = triangle.Vertices[2].Position - boundingBox.Center;

	Vector3F e0 = triangle.Vertices[1].Position - triangle.Vertices[0].Position;
	Vector3F e1 = triangle.Vertices[2].Position - triangle.Vertices[1].Position;
	Vector3F e2 = triangle.Vertices[0].Position - triangle.Vertices[2].Position;

	float p0, p1, p2, min, max, rad, fex, fey, fez;

	fex = std::abs(e0.X);
	fey = std::abs(e0.Y);
	fez = std::abs(e0.Z);
	AXISTEST_X01(e0.Z, e0.Y, fez, fey);
	AXISTEST_Y02(e0.Z, e0.X, fez, fex);
	AXISTEST_Z12(e0.Y, e0.X, fey, fex);

	fex = std::abs(e1.X);
	fey = std::abs(e1.Y);
	fez = std::abs(e1.Z);
	AXISTEST_X01(e1.Z, e1.Y, fez, fey);
	AXISTEST_Y02(e1.Z, e1.X, fez, fex);
	AXISTEST_Z0(e1.Y, e1.X, fey, fex);

	fex = std::abs(e2.X);
	fey = std::abs(e2.Y);
	fez = std::abs(e2.Z);
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

	Vector3F normal = e0.Cross(e1), vMin, vMax;
	for (int i = 0; i < 3; ++i)
	{
		if (normal[i] > 0)
		{
			vMin[i] = -boundingBox.Halfsize[i] - v0[i];
			vMax[i] = boundingBox.Halfsize[i] - v0[i];
		}
		else
		{
			vMin[i] = boundingBox.Halfsize[i] - v0[i];
			vMax[i] = -boundingBox.Halfsize[i] - v0[i];
		}
	}

	if (normal.Dot(vMin) > 0 || normal.Dot(vMax) < 0)
		return false;

	return true;
}