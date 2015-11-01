#pragma once

#include <vector>
#include "Triangle.h"
#include "BoundingBox.h"
#include "Ray.h"

class BVHTreeNode
{
public:
	BVHTreeNode() : left(nullptr), right(nullptr) { };
	~BVHTreeNode() { };
	Triangle* Query(const Ray& ray, float& t) const;

	std::vector<Triangle*> triangles;
	BoundingBox bb;
	BVHTreeNode* left;
	BVHTreeNode* right;
};

class BVHTree
{
public:
	BVHTreeNode* root;

	BVHTree(const std::vector<Triangle*>& triangles, unsigned int minTriangles, unsigned int maxDepth);

	~BVHTree();

	BVHTreeNode* CreateNodeX(const std::vector<Triangle*>& triangles, unsigned int minTriangles, unsigned int maxDepth);
	BVHTreeNode* CreateNodeY(const std::vector<Triangle*>& triangles, unsigned int minTriangles, unsigned int maxDepth);
	BVHTreeNode* CreateNodeZ(const std::vector<Triangle*>& triangles, unsigned int minTriangles, unsigned int maxDepth);
	BVHTreeNode* SplitNode(const std::vector<Triangle*>& triangles, unsigned int minTriangles, unsigned int maxDepth, BoundingBox bb);
	Triangle* Query(const Ray& ray, float& t) const;
};