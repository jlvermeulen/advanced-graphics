#pragma once

#include <vector>
#include "Triangle.h"
#include "BoundingBox.h"
#include "Ray.h"
#include "SpatialIndexDefines.h"

class BVHNode
{
public:
	virtual ~BVHNode() { }
	virtual Triangle* Query(const Ray& ray, float& t) const = 0;
	static BVHNode* Construct(const std::vector<Triangle*>& triangles, const BoundingBox& bb);
	virtual void Compact() { }
	virtual BVHNode** GatherChildren() { return nullptr; }

	BoundingBox bb;
};

class BVHInternal : public BVHNode
{
public:
	BVHInternal() { for (unsigned int i = 0; i < NROFLANES; i++) children[i] = nullptr; }
	~BVHInternal();
	Triangle* Query(const Ray& ray, float& t) const;
	void Compact();
	BVHNode** GetChildren() { return &children[0]; }

	BVHNode* children[NROFLANES];
};

__declspec(align(32))
class BVHLeaf : public BVHNode
{
public:
	BVHLeaf(const std::vector<Triangle*>& triangles, const BoundingBox& bb);
	Triangle* Query(const Ray& ray, float& t) const;

	// crazy shit to ensure alignment
	void* operator new(size_t i){ return _mm_malloc(i, 32); }
	void operator delete(void* p) { _mm_free(p); }

protected:
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
};

class BVHTree
{
public:
	BVHTree();
	BVHTree(const std::vector<Triangle*>& triangles);
	~BVHTree();

	Triangle* Query(const Ray& ray, float& t) const;

private:
	BVHNode* root;

	void Compact();
};