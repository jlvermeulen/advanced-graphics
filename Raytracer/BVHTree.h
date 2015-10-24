#pragma once
#include <deque>
#include "Triangle.h"
#include "BoundingBox.h"
#include "Ray.h"
class BVHTreeNode
{
public:
	BVHTreeNode() : left(nullptr),right(nullptr){};
	~BVHTreeNode(){};
	bool Query(const Ray& ray, Triangle& triangle, double& t) const;

	std::deque<Triangle> triangles;
	BoundingBox bb;
	BVHTreeNode* left;
	BVHTreeNode* right;
};


class BVHTree
{
public:
	BVHTreeNode* root;

	BVHTree(const std::deque<Triangle>& triangles, unsigned int minTriangles, unsigned int maxDepth);

	~BVHTree();

	BVHTreeNode* CreateNodeX(const std::deque<Triangle>& triangles, unsigned int minTriangles, unsigned int maxDepth, BoundingBox bb);
	BVHTreeNode* CreateNodeY(const std::deque<Triangle>& triangles, unsigned int minTriangles, unsigned int maxDepth, BoundingBox bb);
	BVHTreeNode* CreateNodeZ(const std::deque<Triangle>& triangles, unsigned int minTriangles, unsigned int maxDepth, BoundingBox bb);

	bool Query(const Ray& ray, Triangle& triangle, double& t) const;
};

