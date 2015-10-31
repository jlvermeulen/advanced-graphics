#pragma once

#include <vector>
#include "Triangle.h"
#include "BoundingBox.h"
#include "Ray.h"

#define MAXSIZE 128
#define NROFLANES 8

class OctreeNode
{
	friend class GLWidget;

public:
	virtual Triangle* Query(const Ray& ray, float& t) const = 0;

protected:
	BoundingBox bb;
};

class OctreeInternal : public OctreeNode
{
public:
	OctreeInternal(const std::vector<Triangle*>& triangles, const BoundingBox& bb, unsigned int minTriangles, unsigned int maxDepth);
	~OctreeInternal();
	Triangle* Query(const Ray& ray, float& t) const;

	OctreeNode** children;
};

__declspec(align(32))
class OctreeLeaf : public OctreeNode
{
public:
	OctreeLeaf(const std::vector<Triangle*>& triangles, const BoundingBox& bb);
	Triangle* Query(const Ray& ray, float& t) const;

	int count;
	Triangle* triangles[MAXSIZE];

	// Vertices[0]
	union { float vert0X[MAXSIZE]; __m256 vert0X8[MAXSIZE / NROFLANES]; };
	union { float vert0Y[MAXSIZE]; __m256 vert0Y8[MAXSIZE / NROFLANES]; };
	union { float vert0Z[MAXSIZE]; __m256 vert0Z8[MAXSIZE / NROFLANES]; };

	// Vertices[1] - Vertices[0]
	union { float edge1X[MAXSIZE]; __m256 edge1X8[MAXSIZE / NROFLANES]; };
	union { float edge1Y[MAXSIZE]; __m256 edge1Y8[MAXSIZE / NROFLANES]; };
	union { float edge1Z[MAXSIZE]; __m256 edge1Z8[MAXSIZE / NROFLANES]; };

	// Vertices[2] - Vertices[0]
	union { float edge2X[MAXSIZE]; __m256 edge2X8[MAXSIZE / NROFLANES]; };
	union { float edge2Y[MAXSIZE]; __m256 edge2Y8[MAXSIZE / NROFLANES]; };
	union { float edge2Z[MAXSIZE]; __m256 edge2Z8[MAXSIZE / NROFLANES]; };

	// crazy shit to ensure alignment
	void* operator new(size_t i) { return _mm_malloc(i, 32); }
	void operator delete(void* p) { _mm_free(p); }
};

class Octree
{
	friend class GLWidget;

public:
	Octree();
	Octree(const std::vector<Triangle*>& triangles, int minTriangles, int maxDepth);
	~Octree();

	Triangle* Query(const Ray& ray, float& t) const;

private:
	OctreeNode* root;
};