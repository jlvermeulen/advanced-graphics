#include "BVHTree.h"
#include "BoundingBox.h"
#include "Intersections.h"

BVHTree::BVHTree(const std::deque<Triangle>& triangles, unsigned int minTriangles, unsigned int maxDepth)
{
	root = CreateNodeX(triangles, minTriangles, maxDepth, BoundingSphere::FromTriangles(triangles));
}

BVHTree::~BVHTree()
{
}

BVHTreeNode* BVHTree::CreateNodeX(const std::deque<Triangle>& triangles, unsigned int minTriangles, unsigned int maxDepth, BoundingSphere bs)
{
	BVHTreeNode* res = new BVHTreeNode();
	res->bs = bs;
	if (triangles.size() <= minTriangles || maxDepth == 0)
	{
		res->triangles = triangles;
		res->left = nullptr;
		res->right = nullptr;
		return res;
	}
	std::deque<Triangle> leftTriangles, rightTriangles;

	for (Triangle triangle:triangles)
	{
		if (triangle.Vertices[0].Position.X<=bs.center.X)
			leftTriangles.push_back(triangle);
		else
			rightTriangles.push_back(triangle);
	}

	res->left = CreateNodeY(leftTriangles, minTriangles, maxDepth-1, BoundingSphere::FromTriangles(leftTriangles));
	res->right = CreateNodeY(rightTriangles, minTriangles, maxDepth-1, BoundingSphere::FromTriangles(rightTriangles));
	return res;
}

BVHTreeNode* BVHTree::CreateNodeY(const std::deque<Triangle>& triangles, unsigned int minTriangles, unsigned int maxDepth, BoundingSphere bs)
{
	BVHTreeNode* res = new BVHTreeNode();
	res->bs = bs;
	if (triangles.size() <= minTriangles || maxDepth == 0)
	{
		res->triangles = triangles;
		res->left = nullptr;
		res->right = nullptr;
		return res;
	}
	std::deque<Triangle> leftTriangles, rightTriangles;

	for (Triangle triangle : triangles)
	{
		if (triangle.Vertices[0].Position.Y <= bs.center.Y)
			leftTriangles.push_back(triangle);
		else
			rightTriangles.push_back(triangle);
	}

	res->left = CreateNodeZ(leftTriangles, minTriangles, maxDepth - 1, BoundingSphere::FromTriangles(leftTriangles));
	res->right = CreateNodeZ(rightTriangles, minTriangles, maxDepth - 1, BoundingSphere::FromTriangles(rightTriangles));
	return res;
}

BVHTreeNode* BVHTree::CreateNodeZ(const std::deque<Triangle>& triangles, unsigned int minTriangles, unsigned int maxDepth, BoundingSphere bs)
{
	BVHTreeNode* res = new BVHTreeNode();
	res->bs = bs;
	if (triangles.size() <= minTriangles || maxDepth == 0)
	{
		res->triangles = triangles;
		res->left = nullptr;
		res->right = nullptr;
		return res;
	}
	std::deque<Triangle> leftTriangles, rightTriangles;

	for (Triangle triangle : triangles)
	{
		if (triangle.Vertices[0].Position.Z <= bs.center.Z)
			leftTriangles.push_back(triangle);
		else
			rightTriangles.push_back(triangle);
	}

	res->left = CreateNodeX(leftTriangles, minTriangles, maxDepth - 1, BoundingSphere::FromTriangles(leftTriangles));
	res->right = CreateNodeX(rightTriangles, minTriangles, maxDepth - 1, BoundingSphere::FromTriangles(rightTriangles));
	return res;
}


bool BVHTreeNode::Query(const Ray& ray, Triangle& triangle, double& t) const
{
	double tBox = std::numeric_limits<double>::min();
	if (!Intersects(ray, bs, tBox) || tBox > t)
		return false;

	bool hit = false;
	for (std::deque<Triangle>::const_iterator it = triangles.begin(); it != triangles.end(); ++it)
	{
		double tCurr = 0;
		if (Intersects(ray, *it, tCurr) && tCurr < t)
		{
			t = tCurr;
			triangle = *it;
			hit = true;
		}
	}

	if (left != nullptr)
		hit |= left->Query(ray, triangle, t);
	if (right != nullptr)
		hit |= right->Query(ray, triangle, t);

	return hit;
}

bool BVHTree::Query(const Ray& ray, Triangle& triangle, double& t) const
{
	t = std::numeric_limits<double>::max();
	return root->Query(ray, triangle, t);
}