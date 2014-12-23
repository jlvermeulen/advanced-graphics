#pragma once

#include <vector>
#include <deque>
#include "Triangle.h"
#include "BoundingBox.h"
#include "Ray.h"

class OctreeNode
{
public:
	OctreeNode();
	OctreeNode(const std::deque<Triangle>& triangles, const BoundingBox& bb, unsigned int minTriangles, unsigned int maxDepth);

	bool Query(const Ray& ray, Triangle& triangle, double& t) const;

	BoundingBox bb;
private:
	std::deque<Triangle> triangles;
	std::vector<OctreeNode> children;
};

class Octree
{
public:
	Octree();
	Octree(const std::deque<Triangle>& triangles, int minTriangles, int maxDepth);

	bool Query(const Ray& ray, Triangle& triangle, double& t) const;

	OctreeNode root;
private:
};