#include "Octree.h"
#include "Intersections.h"
#include <limits>
#include <immintrin.h>

Octree::Octree() { }
Octree::Octree(const std::vector<Triangle*>& triangles, int minTriangles, int maxDepth)
{
	if (triangles.size() <= MAXSIZE || maxDepth == 0)
		root = new OctreeLeaf(triangles, BoundingBox::FromTriangles(triangles));
	else
		root = new OctreeInternal(triangles, BoundingBox::FromTriangles(triangles), minTriangles, maxDepth);
}
Octree::~Octree() { delete root; }

Triangle* Octree::Query(const Ray& ray, double& t) const
{
	t = std::numeric_limits<double>::max();
	return root->Query(ray, t);
}

OctreeInternal::OctreeInternal(const std::vector<Triangle*>& triangles, const BoundingBox& bb, unsigned int minTriangles, unsigned int maxDepth)
{
	this->bb = bb;
	Vector3D half = bb.Halfsize / 2;
	children = new OctreeNode*[8];
	for (int i = 0; i < 8; ++i)
	{
		double x = bb.Center.X + (i & 4 ? half.X : -half.X);
		double y = bb.Center.Y + (i & 2 ? half.Y : -half.Y);
		double z = bb.Center.Z + (i & 1 ? half.Z : -half.Z);
		BoundingBox childBB(Vector3D(x, y, z), half);

		std::vector<Triangle*> childTriangles;
		for (int j = 0; j < triangles.size(); j++)
			if (Intersects(*triangles[j], childBB))
				childTriangles.push_back(triangles[j]);

		if (childTriangles.size() <= MAXSIZE || maxDepth == 0)
			children[i] = new OctreeLeaf(childTriangles, childBB);
		else
			children[i] = new OctreeInternal(childTriangles, childBB, minTriangles, maxDepth - 1);
	}
}

OctreeInternal::~OctreeInternal()
{
	if (children == nullptr)
		return;

	for (int i = 0; i < 8; ++i)
		delete children[i];
	delete [] children;
}

Triangle* OctreeInternal::Query(const Ray& ray, double& t) const
{
	double tBox = std::numeric_limits<double>::min();
	if (!Intersects(ray, bb, tBox) || tBox > t)
		return nullptr;

	Triangle* triangle = nullptr;
	for (int i = 0; i < 8; ++i)
	{
		Triangle* tri = children[i]->Query(ray, t);
		triangle = tri == nullptr ? triangle : tri;
	}
	return triangle;
}

OctreeLeaf::OctreeLeaf(const std::vector<Triangle*>& triangles, const BoundingBox& bb)
{
	this->bb = bb;

	int i = 0;
	for (; i < triangles.size(); i++)
	{
		this->triangles[i] = triangles[i];

		this->vert0X[i] = triangles[i]->Vertices[0].Position.X;
		this->vert0Y[i] = triangles[i]->Vertices[0].Position.Y;
		this->vert0Z[i] = triangles[i]->Vertices[0].Position.Z;

		this->edge1X[i] = triangles[i]->Vertices[1].Position.X - triangles[i]->Vertices[0].Position.X;
		this->edge1Y[i] = triangles[i]->Vertices[1].Position.Y - triangles[i]->Vertices[0].Position.Y;
		this->edge1Z[i] = triangles[i]->Vertices[1].Position.Z - triangles[i]->Vertices[0].Position.Z;

		this->edge2X[i] = triangles[i]->Vertices[2].Position.X - triangles[i]->Vertices[0].Position.X;
		this->edge2Y[i] = triangles[i]->Vertices[2].Position.Y - triangles[i]->Vertices[0].Position.Y;
		this->edge2Z[i] = triangles[i]->Vertices[2].Position.Z - triangles[i]->Vertices[0].Position.Z;
	}

	for (; i < MAXSIZE; i++)
	{
		vert0X[i] = vert0Y[i] = vert0Z[i] = 0;
		edge1X[i] = edge1Y[i] = edge1Z[i] = 0;
		edge2X[i] = edge2Y[i] = edge2Z[i] = 0;
	}

	this->count = triangles.size() / 4;
	if (triangles.size() % 4 != 0)
		this->count++;
}

Triangle* OctreeLeaf::Query(const Ray& ray, double& t) const
{
	double tBox = std::numeric_limits<double>::min();
	if (!Intersects(ray, bb, tBox) || tBox > t)
		return nullptr;

	const __m256d rayDirX = _mm256_set1_pd(ray.Direction.X);
	const __m256d rayDirY = _mm256_set1_pd(ray.Direction.Y);
	const __m256d rayDirZ = _mm256_set1_pd(ray.Direction.Z);

	const __m256d rayPosX = _mm256_set1_pd(ray.Origin.X);
	const __m256d rayPosY = _mm256_set1_pd(ray.Origin.Y);
	const __m256d rayPosZ = _mm256_set1_pd(ray.Origin.Z);

	union { __m256d distances[MAXSIZE / NROFLANES]; double dists[MAXSIZE]; };

	for (int i = 0; i < count; i++)
	{
		// Vector3D e1 = triangle.Vertices[1].Position - triangle.Vertices[0].Position;
		const __m256d e1X = edge1X4[i];
		const __m256d e1Y = edge1Y4[i];
		const __m256d e1Z = edge1Z4[i];

		// Vector3D e2 = triangle.Vertices[2].Position - triangle.Vertices[0].Position;
		const __m256d e2X = edge2X4[i];
		const __m256d e2Y = edge2Y4[i];
		const __m256d e2Z = edge2Z4[i];

		// Vector3D p = ray.Direction.Cross(e2);
		const __m256d pX = _mm256_sub_pd(_mm256_mul_pd(rayDirY, e2Z), _mm256_mul_pd(rayDirZ, e2Y));
		const __m256d pY = _mm256_sub_pd(_mm256_mul_pd(rayDirZ, e2X), _mm256_mul_pd(rayDirX, e2Z));
		const __m256d pZ = _mm256_sub_pd(_mm256_mul_pd(rayDirX, e2Y), _mm256_mul_pd(rayDirY, e2X));

		// double det = e1.Dot(p);
		const __m256d det = _mm256_add_pd(_mm256_mul_pd(e1X, pX), _mm256_add_pd(_mm256_mul_pd(e1Y, pY), _mm256_mul_pd(e1Z, pZ)));

		// if (det > -EPSILON && det < EPSILON)
		//     return false;
		__m256d mask = _mm256_or_pd(_mm256_cmp_pd(det, _mm256_set1_pd(-EPSILON), _CMP_LE_OS), _mm256_cmp_pd(det, _mm256_set1_pd(EPSILON), _CMP_GE_OS));

		// double invDet = 1 / det;
		const __m256d invDet = _mm256_div_pd(_mm256_set1_pd(1.0), det);

		// Vector3D r = ray.Origin - triangle.Vertices[0].Position;
		const __m256d rX = _mm256_sub_pd(rayPosX, vert0X4[i]);
		const __m256d rY = _mm256_sub_pd(rayPosY, vert0Y4[i]);
		const __m256d rZ = _mm256_sub_pd(rayPosZ, vert0Z4[i]);

		// double u = r.Dot(p) * invDet;
		const __m256d u = _mm256_mul_pd(invDet, _mm256_add_pd(_mm256_mul_pd(rX, pX), _mm256_add_pd(_mm256_mul_pd(rY, pY), _mm256_mul_pd(rZ, pZ))));

		// if (u < 0 || u > 1)
		//	   return false;
		mask = _mm256_and_pd(mask, _mm256_cmp_pd(u, _mm256_setzero_pd(), _CMP_GE_OS));

		// Vector3D q = r.Cross(e1);
		const __m256d qX = _mm256_sub_pd(_mm256_mul_pd(rY, e1Z), _mm256_mul_pd(rZ, e1Y));
		const __m256d qY = _mm256_sub_pd(_mm256_mul_pd(rZ, e1X), _mm256_mul_pd(rX, e1Z));
		const __m256d qZ = _mm256_sub_pd(_mm256_mul_pd(rX, e1Y), _mm256_mul_pd(rY, e1X));

		// double v = ray.Direction.Dot(q) * invDet;
		const __m256d v = _mm256_mul_pd(invDet, _mm256_add_pd(_mm256_mul_pd(rayDirX, qX), _mm256_add_pd(_mm256_mul_pd(rayDirY, qY), _mm256_mul_pd(rayDirZ, qZ))));

		// if (v < 0 || u + v > 1)
		//     return false;
		mask = _mm256_and_pd(mask, _mm256_and_pd(_mm256_cmp_pd(v, _mm256_setzero_pd(), _CMP_GE_OS), _mm256_cmp_pd(_mm256_add_pd(u, v), _mm256_set1_pd(1.0), _CMP_LE_OS)));

		// double tt = e2.Dot(q) * invDet;
		const __m256d tt = _mm256_mul_pd(invDet, _mm256_add_pd(_mm256_mul_pd(e2X, qX), _mm256_add_pd(_mm256_mul_pd(e2Y, qY), _mm256_mul_pd(e2Z, qZ))));

		// if (tt > EPSILON)
		// {
		//     t = tt;
		//     return true;
		// }
		//
		// return false;
		distances[i] = _mm256_and_pd(tt, mask);
	}

	Triangle* triangle = nullptr;
	for (int i = 0; i < count * 4; i++)
		if (dists[i] < t && dists[i] > EPSILON)
		{
			t = dists[i];
			triangle = triangles[i];
		}

	return triangle;
}