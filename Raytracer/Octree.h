#pragma once

#include <vector>
#include "Triangle.h"
#include "BoundingBox.h"
#include "Ray.h"

class OctreeNode
{
public:
	OctreeNode();
	OctreeNode(const std::vector<Triangle>& triangles, const BoundingBox& bb, unsigned int minTriangles, unsigned int maxDepth);
	~OctreeNode();

	bool Query(const Ray& ray, Triangle& triangle, double& t) const;

private:
	std::vector<Triangle> triangles;
	OctreeNode* children;
	BoundingBox bb;
};

class Octree
{
public:
	Octree();
	Octree(const std::vector<Triangle>& triangles, int minTriangles, int maxDepth);

	bool Query(const Ray& ray, Triangle& triangle, double& t) const;

private:
	OctreeNode root;
};