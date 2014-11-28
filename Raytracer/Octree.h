#pragma once

#include <vector>
#include "TriangleD.h"
#include "BoundingBox.h"

class OctreeNode
{
public:
	OctreeNode();
	OctreeNode(const std::vector<TriangleD>& triangles, const BoundingBox& bb, unsigned int minTriangles, unsigned int maxDepth);
	~OctreeNode();

private:
	std::vector<TriangleD> triangles;
	OctreeNode* children;
	BoundingBox bb;
};

class Octree
{
public:
	Octree();
	Octree(const std::vector<TriangleD>& triangles, int minTriangles, int maxDepth);

private:
	OctreeNode root;
};