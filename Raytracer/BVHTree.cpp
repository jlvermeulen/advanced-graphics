#include "BVHTree.h"
#include "BoundingBox.h"
#include "Intersections.h"
#include <stdexcept>

BVHTree::BVHTree(const std::deque<Triangle>& triangles, unsigned int minTriangles, unsigned int maxDepth)
{
	root = SplitNode(triangles, minTriangles, maxDepth, BoundingBox::FromTriangles(triangles));
}

BVHTree::~BVHTree()
{
}

BVHTreeNode* BVHTree::CreateNodeX(const std::deque<Triangle>& triangles, unsigned int minTriangles, unsigned int maxDepth)
{
	BoundingBox bb = BoundingBox::FromTriangles(triangles);
	BVHTreeNode* res = new BVHTreeNode();
	res->bb = bb;
	if (triangles.size() <= minTriangles || maxDepth == 0)
	{
		res->triangles = triangles;
		res->left = nullptr;
		res->right = nullptr;
		return res;
	}

	Vector3D one18X = Vector3D(bb.Halfsize.X/9, 0, 0); // used to find other halfsizes + centers

	BoundingBox bbl[8];// left halfs
	BoundingBox bbr[8];//  right halfs

	// create partitionings
	for (int i = 0; i < 8; ++i)
	{
		bbl[i] = BoundingBox(bb.Center - (8 - i) * one18X, bb.Halfsize - (8 - i)*one18X);
		bbr[i] = BoundingBox(bb.Center + (i+1)*one18X, bb.Halfsize - (1 + i)*one18X);
	}

	std::deque<Triangle> leftTriangles[8], rightTriangles[8];

	// add triangles to partitions
	for (Triangle triangle : triangles)
	{
		for (int i = 0; i < 8; ++i)
		{
			if (Intersects(triangle, bbl[i]))
				leftTriangles[i].push_back(triangle);
			else
				rightTriangles[i].push_back(triangle);
		}
	}


	// choose best one
	double costNoSplit = bb.SurfaceArea()*triangles.size();
	std::pair<double, int> best(costNoSplit, -1);
	for (int i = 0; i < 8; ++i)
	{
		bbl[i] = BoundingBox::FromTriangles(leftTriangles[i]);
		bbr[i] = BoundingBox::FromTriangles(rightTriangles[i]);
		double score = bbl[i].SurfaceArea() * leftTriangles[i].size() + bbr[i].SurfaceArea() * rightTriangles[i].size();
		if (score < best.first)
			best = std::pair<double, int>(score, i);
	}
	if (best.second == -1)
	{
		res->triangles = triangles;
		res->left = nullptr;
		res->right = nullptr;
		return res;
	}

	res->left = SplitNode(leftTriangles[best.second], minTriangles, maxDepth - 1, bbl[best.second]);
	res->right = SplitNode(rightTriangles[best.second], minTriangles, maxDepth - 1, bbr[best.second]);
	return res;
}

BVHTreeNode* BVHTree::CreateNodeY(const std::deque<Triangle>& triangles, unsigned int minTriangles, unsigned int maxDepth)
{
	BoundingBox bb = BoundingBox::FromTriangles(triangles);
	BVHTreeNode* res = new BVHTreeNode();
	res->bb = bb;
	if (triangles.size() <= minTriangles || maxDepth == 0)
	{
		res->triangles = triangles;
		res->left = nullptr;
		res->right = nullptr;
		return res;
	}

	Vector3D one18Y = Vector3D(0,bb.Halfsize.Y/9, 0); // used to find other halfsizes + centers

	BoundingBox bbl[8];// left halfs
	BoundingBox bbr[8];//  right halfs

	// create partitionings
	for (int i = 0; i < 8; ++i)
	{
		bbl[i] = BoundingBox(bb.Center - (8 - i) * one18Y, bb.Halfsize - (8 - i)*one18Y);
		bbr[i] = BoundingBox(bb.Center + (i+1)*one18Y, bb.Halfsize - (1 + i)*one18Y);
	}

	std::deque<Triangle> leftTriangles[8], rightTriangles[8];

	// add triangles to partitions
	for (Triangle triangle : triangles)
	{
		for (int i = 0; i < 8; ++i)
		{
			if (Intersects(triangle, bbl[i]))
				leftTriangles[i].push_back(triangle);
			else
				rightTriangles[i].push_back(triangle);
		}
	}

	// choose best one
	double costNoSplit = bb.SurfaceArea()*triangles.size();
	std::pair<double, int> best(costNoSplit, -1);
	for (int i = 0; i < 8; ++i)
	{
		bbl[i] = BoundingBox::FromTriangles(leftTriangles[i]);
		bbr[i] = BoundingBox::FromTriangles(rightTriangles[i]);
		double score = bbl[i].SurfaceArea() * leftTriangles[i].size() + bbr[i].SurfaceArea() * rightTriangles[i].size();
		if (score < best.first)
			best = std::pair<double, int>(score, i);
	}
	if (best.second == -1)
	{
		res->triangles = triangles;
		res->left = nullptr;
		res->right = nullptr;
		return res;
	}

	res->left = SplitNode(leftTriangles[best.second], minTriangles, maxDepth - 1, bbl[best.second]);
	res->right = SplitNode(rightTriangles[best.second], minTriangles, maxDepth - 1, bbr[best.second]);
	return res;
}

BVHTreeNode* BVHTree::CreateNodeZ(const std::deque<Triangle>& triangles, unsigned int minTriangles, unsigned int maxDepth)
{
	BoundingBox bb = BoundingBox::FromTriangles(triangles);
	BVHTreeNode* res = new BVHTreeNode();
	res->bb = bb;
	if (triangles.size() <= minTriangles || maxDepth == 0)
	{
		res->triangles = triangles;
		res->left = nullptr;
		res->right = nullptr;
		return res;
	}

	Vector3D one18Z = Vector3D(0,0,bb.Halfsize.Z/9); // used to find other halfsizes + centers

	BoundingBox bbl[8];// left halfs
	BoundingBox bbr[8];//  right halfs

	// create partitionings
	for (int i = 0; i < 8; ++i)
	{
		bbl[i] = BoundingBox(bb.Center - (8 - i) * one18Z, bb.Halfsize - (8 - i)*one18Z);
		bbr[i] = BoundingBox(bb.Center + (1+i)*one18Z, bb.Halfsize - (1 + i)*one18Z);
	}

	std::deque<Triangle> leftTriangles[8], rightTriangles[8];

	// add triangles to partitions
	for (Triangle triangle : triangles)
	{
		for (int i = 0; i < 8; ++i)
		{
			if (Intersects(triangle, bbl[i]))
				leftTriangles[i].push_back(triangle);
			else
				rightTriangles[i].push_back(triangle);
		}
	}

	// choose best one
	double costNoSplit = bb.SurfaceArea()*triangles.size();
	std::pair<double, int> best(costNoSplit, -1);
	for (int i = 0; i < 8; ++i)
	{
		bbl[i] = BoundingBox::FromTriangles(leftTriangles[i]);
		bbr[i] = BoundingBox::FromTriangles(rightTriangles[i]);
		double score = bbl[i].SurfaceArea() * leftTriangles[i].size() + bbr[i].SurfaceArea() * rightTriangles[i].size();
		if (score < best.first)
			best = std::pair<double, int>(score, i);
	}
	if (best.second == -1)
	{
		res->triangles = triangles;
		res->left = nullptr;
		res->right = nullptr;
		return res;
	}

	res->left = SplitNode(leftTriangles[best.second], minTriangles, maxDepth - 1, bbl[best.second]);
	res->right = SplitNode(rightTriangles[best.second], minTriangles, maxDepth - 1, bbr[best.second]);
	return res;
}


BVHTreeNode* BVHTree::SplitNode(const std::deque<Triangle>& triangles, unsigned int minTriangles, unsigned int maxDepth, BoundingBox bb)
{
	if (bb.Halfsize.X >= bb.Halfsize.Y)
		if (bb.Halfsize.X > bb.Halfsize.Z)
			return CreateNodeX(triangles, minTriangles, maxDepth);

	if (bb.Halfsize.Y >= bb.Halfsize.Z)
		return CreateNodeY(triangles, minTriangles, maxDepth);
	else
		return CreateNodeZ(triangles, minTriangles, maxDepth);

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