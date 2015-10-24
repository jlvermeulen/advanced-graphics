#include "BVHTree.h"
#include "BoundingBox.h"
#include "Intersections.h"

BVHTree::BVHTree(const std::deque<Triangle>& triangles, unsigned int minTriangles, unsigned int maxDepth)
{
	root = CreateNodeX(triangles, minTriangles, maxDepth, BoundingBox::FromTriangles(triangles));
}

BVHTree::~BVHTree()
{
}

BVHTreeNode* BVHTree::CreateNodeX(const std::deque<Triangle>& triangles, unsigned int minTriangles, unsigned int maxDepth, BoundingBox bb)
{
	BVHTreeNode* res = new BVHTreeNode();
		res->bb = bb;
	if (triangles.size() <= minTriangles || maxDepth == 0)
	{
		res->triangles = triangles;
		res->left = nullptr;
		res->right = nullptr;
		return res;
	}
	Vector3D leftCenter = bb.Center,rightCenter = bb.Center;
	Vector3D quarter = bb.Halfsize;
	quarter.X /= 2;
	leftCenter.X -= quarter.X;
	rightCenter.X += quarter.X;

	BoundingBox leftBB = BoundingBox(leftCenter, quarter), rightBB = BoundingBox(rightCenter, quarter);
	std::deque<Triangle> leftTriangles, rightTriangles;

	for (Triangle triangle:triangles)
	{
		if (Intersects(triangle, leftBB))
			leftTriangles.push_back(triangle);
		if (Intersects(triangle, rightBB))
			rightTriangles.push_back(triangle);
	}

	res->left = CreateNodeY(leftTriangles, minTriangles, maxDepth-1, leftBB);
	res->right = CreateNodeY(rightTriangles, minTriangles, maxDepth-1, rightBB);
	return res;
}

BVHTreeNode* BVHTree::CreateNodeY(const std::deque<Triangle>& triangles, unsigned int minTriangles, unsigned int maxDepth, BoundingBox bb)
{
	BVHTreeNode* res = new BVHTreeNode();
	res->bb = bb;
	if (triangles.size() <= minTriangles || maxDepth == 0)
	{
		res->triangles = triangles;
		res->left = nullptr;
		res->right = nullptr;
		return res;
	}
	Vector3D leftCenter = bb.Center, rightCenter = bb.Center;
	Vector3D quarter = bb.Halfsize;
	quarter.Y /= 2;
	leftCenter.Y -= quarter.Y;
	rightCenter.Y += quarter.Y;

	BoundingBox leftBB = BoundingBox(leftCenter, quarter), rightBB = BoundingBox(rightCenter, quarter);
	std::deque<Triangle> leftTriangles, rightTriangles;

	for (Triangle triangle : triangles)
	{
		if (Intersects(triangle, leftBB))
			leftTriangles.push_back(triangle);
		if (Intersects(triangle, rightBB))
			rightTriangles.push_back(triangle);
	}

	res->left = CreateNodeZ(leftTriangles, minTriangles, maxDepth - 1, leftBB);
	res->right = CreateNodeZ(rightTriangles, minTriangles, maxDepth - 1, rightBB);
	return res;
}

BVHTreeNode* BVHTree::CreateNodeZ(const std::deque<Triangle>& triangles, unsigned int minTriangles, unsigned int maxDepth, BoundingBox bb)
{
	BVHTreeNode* res = new BVHTreeNode();
	res->bb = bb;
	if (triangles.size() <= minTriangles || maxDepth == 0)
	{
		res->triangles = triangles;
		res->left = nullptr;
		res->right = nullptr;
		return res;
	}
	Vector3D leftCenter = bb.Center, rightCenter = bb.Center;
	Vector3D quarter = bb.Halfsize;
	quarter.Z /= 2;
	leftCenter.Z -= quarter.Z;
	rightCenter.Z += quarter.Z;

	BoundingBox leftBB = BoundingBox(leftCenter, quarter), rightBB = BoundingBox(rightCenter, quarter);
	std::deque<Triangle> leftTriangles, rightTriangles;

	for (Triangle triangle : triangles)
	{
		if (Intersects(triangle, leftBB))
			leftTriangles.push_back(triangle);
		if (Intersects(triangle, rightBB))
			rightTriangles.push_back(triangle);
	}

	res->left = CreateNodeX(leftTriangles, minTriangles, maxDepth - 1, leftBB);
	res->right = CreateNodeX(rightTriangles, minTriangles, maxDepth - 1, rightBB);
	return res;
}


bool BVHTreeNode::Query(const Ray& ray, Triangle& triangle, double& t) const
{
	double tBox = std::numeric_limits<double>::min();
	if (!Intersects(ray, bb, tBox) || tBox > t)
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