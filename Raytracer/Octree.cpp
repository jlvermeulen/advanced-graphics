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

Triangle* Octree::Query(const Ray& ray, float& t) const
{
	t = std::numeric_limits<float>::max();
	return root->Query(ray, t);
}

OctreeInternal::OctreeInternal(const std::vector<Triangle*>& triangles, const BoundingBox& bb, unsigned int minTriangles, unsigned int maxDepth)
{
	this->bb = bb;
	Vector3F half = bb.Halfsize / 2;
	children = new OctreeNode*[8];
	for (int i = 0; i < 8; ++i)
	{
		float x = bb.Center.X + (i & 4 ? half.X : -half.X);
		float y = bb.Center.Y + (i & 2 ? half.Y : -half.Y);
		float z = bb.Center.Z + (i & 1 ? half.Z : -half.Z);
		BoundingBox childBB(Vector3F(x, y, z), half);

		std::vector<Triangle*> childTriangles;
		for (unsigned int j = 0; j < triangles.size(); j++)
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

Triangle* OctreeInternal::Query(const Ray& ray, float& t) const
{
	float tBox = std::numeric_limits<float>::min();
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

	unsigned int i = 0;
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
		vert0X[i] = vert0Y[i] = vert0Z[i] = 0.0f;
		edge1X[i] = edge1Y[i] = edge1Z[i] = 0.0f;
		edge2X[i] = edge2Y[i] = edge2Z[i] = 0.0f;
	}

	this->count = triangles.size() / NROFLANES;
	if (triangles.size() % NROFLANES != 0)
		this->count++;
}

Triangle* OctreeLeaf::Query(const Ray& ray, float& t) const
{
	float tBox = std::numeric_limits<float>::min();
	if (!Intersects(ray, bb, tBox) || tBox > t)
		return nullptr;

	const __m256 rayDirX = _mm256_set1_ps(ray.Direction.X);
	const __m256 rayDirY = _mm256_set1_ps(ray.Direction.Y);
	const __m256 rayDirZ = _mm256_set1_ps(ray.Direction.Z);

	const __m256 rayPosX = _mm256_set1_ps(ray.Origin.X);
	const __m256 rayPosY = _mm256_set1_ps(ray.Origin.Y);
	const __m256 rayPosZ = _mm256_set1_ps(ray.Origin.Z);

	union { float dists[MAXSIZE]; __m256 distances[MAXSIZE / NROFLANES]; };

	for (int i = 0; i < count; i++)
	{
		// Vector3F e1 = triangle.Vertices[1].Position - triangle.Vertices[0].Position;
		const __m256 e1X = edge1X8[i];
		const __m256 e1Y = edge1Y8[i];
		const __m256 e1Z = edge1Z8[i];

		// Vector3F e2 = triangle.Vertices[2].Position - triangle.Vertices[0].Position;
		const __m256 e2X = edge2X8[i];
		const __m256 e2Y = edge2Y8[i];
		const __m256 e2Z = edge2Z8[i];

		// Vector3F p = ray.Direction.Cross(e2);
		const __m256 pX = _mm256_sub_ps(_mm256_mul_ps(rayDirY, e2Z), _mm256_mul_ps(rayDirZ, e2Y));
		const __m256 pY = _mm256_sub_ps(_mm256_mul_ps(rayDirZ, e2X), _mm256_mul_ps(rayDirX, e2Z));
		const __m256 pZ = _mm256_sub_ps(_mm256_mul_ps(rayDirX, e2Y), _mm256_mul_ps(rayDirY, e2X));

		// float det = e1.Dot(p);
		const __m256 det = _mm256_add_ps(_mm256_mul_ps(e1X, pX), _mm256_add_ps(_mm256_mul_ps(e1Y, pY), _mm256_mul_ps(e1Z, pZ)));

		// if (det > -EPSILON && det < EPSILON)
		//     return false;
		__m256 mask = _mm256_or_ps(_mm256_cmp_ps(det, _mm256_set1_ps(-EPSILON), _CMP_LE_OS), _mm256_cmp_ps(det, _mm256_set1_ps(EPSILON), _CMP_GE_OS));

		// float invDet = 1 / det;
		const __m256 invDet = _mm256_div_ps(_mm256_set1_ps(1.0f), det);

		// Vector3F r = ray.Origin - triangle.Vertices[0].Position;
		const __m256 rX = _mm256_sub_ps(rayPosX, vert0X8[i]);
		const __m256 rY = _mm256_sub_ps(rayPosY, vert0Y8[i]);
		const __m256 rZ = _mm256_sub_ps(rayPosZ, vert0Z8[i]);

		// float u = r.Dot(p) * invDet;
		const __m256 u = _mm256_mul_ps(invDet, _mm256_add_ps(_mm256_mul_ps(rX, pX), _mm256_add_ps(_mm256_mul_ps(rY, pY), _mm256_mul_ps(rZ, pZ))));

		// if (u < 0 || u > 1)
		//	   return false;
		mask = _mm256_and_ps(mask, _mm256_cmp_ps(u, _mm256_setzero_ps(), _CMP_GE_OS));

		// Vector3F q = r.Cross(e1);
		const __m256 qX = _mm256_sub_ps(_mm256_mul_ps(rY, e1Z), _mm256_mul_ps(rZ, e1Y));
		const __m256 qY = _mm256_sub_ps(_mm256_mul_ps(rZ, e1X), _mm256_mul_ps(rX, e1Z));
		const __m256 qZ = _mm256_sub_ps(_mm256_mul_ps(rX, e1Y), _mm256_mul_ps(rY, e1X));

		// float v = ray.Direction.Dot(q) * invDet;
		const __m256 v = _mm256_mul_ps(invDet, _mm256_add_ps(_mm256_mul_ps(rayDirX, qX), _mm256_add_ps(_mm256_mul_ps(rayDirY, qY), _mm256_mul_ps(rayDirZ, qZ))));

		// if (v < 0 || u + v > 1)
		//     return false;
		mask = _mm256_and_ps(mask, _mm256_and_ps(_mm256_cmp_ps(v, _mm256_setzero_ps(), _CMP_GE_OS), _mm256_cmp_ps(_mm256_add_ps(u, v), _mm256_set1_ps(1.0f), _CMP_LE_OS)));

		// float tt = e2.Dot(q) * invDet;
		const __m256 tt = _mm256_mul_ps(invDet, _mm256_add_ps(_mm256_mul_ps(e2X, qX), _mm256_add_ps(_mm256_mul_ps(e2Y, qY), _mm256_mul_ps(e2Z, qZ))));

		// if (tt > EPSILON)
		// {
		//     t = tt;
		//     return true;
		// }
		//
		// return false;
		distances[i] = _mm256_and_ps(tt, mask);
	}

	Triangle* triangle = nullptr;
	for (int i = 0; i < count * NROFLANES; i++)
		if (dists[i] < t && dists[i] > EPSILON)
		{
			t = dists[i];
			triangle = triangles[i];
		}

	return triangle;
}