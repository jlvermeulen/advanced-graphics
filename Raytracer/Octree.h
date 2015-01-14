#pragma once

#include <deque>
#include "Triangle.h"
#include "BoundingBox.h"
#include "Ray.h"

class OctreeNode
{
	friend class GLWidget;

public:
	OctreeNode();
	OctreeNode(const std::deque<Triangle>& triangles, const BoundingBox& bb, unsigned int minTriangles, unsigned int maxDepth);
	~OctreeNode();

	bool Query(const Ray& ray, Triangle& triangle, double& t) const;

private:
	std::deque<Triangle> triangles;
	OctreeNode** children;
	BoundingBox bb;
};

class Octree
{
	friend class GLWidget;

public:
	Octree();
	Octree(const std::deque<Triangle>& triangles, int minTriangles, int maxDepth);

	bool Query(const Ray& ray, Triangle& triangle, double& t) const;

private:
	OctreeNode root;
};