#pragma once

#include <vector>
#include "Triangle.h"
#include "BoundingBox.h"
#include "Ray.h"

#define MAXSIZE 64
#define NROFLANES 4

class OctreeNode
{
	friend class GLWidget;

public:
	virtual Triangle* Query(const Ray& ray, double& t) const = 0;

protected:
	BoundingBox bb;
};

class OctreeInternal : public OctreeNode
{
public:
	OctreeInternal(const std::vector<Triangle*>& triangles, const BoundingBox& bb, unsigned int minTriangles, unsigned int maxDepth);
	~OctreeInternal();
	Triangle* Query(const Ray& ray, double& t) const;

	OctreeNode** children;
};

class OctreeLeaf : public OctreeNode
{
public:
	OctreeLeaf(const std::vector<Triangle*>& triangles, const BoundingBox& bb);
	Triangle* Query(const Ray& ray, double& t) const;

	int count;
	Triangle* triangles[MAXSIZE];

	// Vertices[0]
	union { double vert0X[MAXSIZE]; __m256d vert0X4[MAXSIZE / NROFLANES]; };
	union { double vert0Y[MAXSIZE]; __m256d vert0Y4[MAXSIZE / NROFLANES]; };
	union { double vert0Z[MAXSIZE]; __m256d vert0Z4[MAXSIZE / NROFLANES]; };

	// Vertices[1] - Vertices[0]
	union { double edge1X[MAXSIZE]; __m256d edge1X4[MAXSIZE / NROFLANES]; };
	union { double edge1Y[MAXSIZE]; __m256d edge1Y4[MAXSIZE / NROFLANES]; };
	union { double edge1Z[MAXSIZE]; __m256d edge1Z4[MAXSIZE / NROFLANES]; };

	// Vertices[2] - Vertices[0]
	union { double edge2X[MAXSIZE]; __m256d edge2X4[MAXSIZE / NROFLANES]; };
	union { double edge2Y[MAXSIZE]; __m256d edge2Y4[MAXSIZE / NROFLANES]; };
	union { double edge2Z[MAXSIZE]; __m256d edge2Z4[MAXSIZE / NROFLANES]; };
};

class Octree
{
	friend class GLWidget;

public:
	Octree();
	Octree(const std::vector<Triangle*>& triangles, int minTriangles, int maxDepth);
	~Octree();

	Triangle* Query(const Ray& ray, double& t) const;

private:
	OctreeNode* root;
};