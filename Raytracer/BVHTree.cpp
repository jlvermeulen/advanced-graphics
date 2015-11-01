#include "BVHTree.h"
#include "Intersections.h"

BVHTree::BVHTree(const std::vector<Triangle*>& triangles, unsigned int minTriangles, unsigned int maxDepth)
{
	root = SplitNode(triangles, minTriangles, maxDepth, BoundingBox::FromTriangles(triangles));
}

BVHTree::~BVHTree()
{
}

BVHTreeNode* BVHTree::CreateNodeX(const std::vector<Triangle*>& triangles, unsigned int minTriangles, unsigned int maxDepth)
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

	Vector3F one18X = Vector3F(bb.Halfsize.X/9, 0, 0); // used to find other halfsizes + centers

	BoundingBox bbl[8];// left halfs
	BoundingBox bbr[8];//  right halfs

	// create partitionings
	for (int i = 0; i < 8; ++i)
	{
		bbl[i] = BoundingBox(bb.Center - (8 - i) * one18X, bb.Halfsize - (8 - i) * one18X);
		bbr[i] = BoundingBox(bb.Center + (i + 1) * one18X, bb.Halfsize - (i + 1) * one18X);
	}

	std::vector<Triangle*> leftTriangles[8], rightTriangles[8];

	// add triangles to partitions
	for (Triangle* triangle : triangles)
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

BVHTreeNode* BVHTree::CreateNodeY(const std::vector<Triangle*>& triangles, unsigned int minTriangles, unsigned int maxDepth)
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

	Vector3F one18Y = Vector3F(0,bb.Halfsize.Y/9, 0); // used to find other halfsizes + centers

	BoundingBox bbl[8];// left halfs
	BoundingBox bbr[8];//  right halfs

	// create partitionings
	for (int i = 0; i < 8; ++i)
	{
		bbl[i] = BoundingBox(bb.Center - (8 - i) * one18Y, bb.Halfsize - (8 - i) * one18Y);
		bbr[i] = BoundingBox(bb.Center + (i + 1) * one18Y, bb.Halfsize - (i + 1) * one18Y);
	}

	std::vector<Triangle*> leftTriangles[8], rightTriangles[8];

	// add triangles to partitions
	for (Triangle* triangle : triangles)
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

BVHTreeNode* BVHTree::CreateNodeZ(const std::vector<Triangle*>& triangles, unsigned int minTriangles, unsigned int maxDepth)
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

	Vector3F one18Z = Vector3F(0,0,bb.Halfsize.Z/9); // used to find other halfsizes + centers

	BoundingBox bbl[8];// left halfs
	BoundingBox bbr[8];//  right halfs

	// create partitionings
	for (int i = 0; i < 8; ++i)
	{
		bbl[i] = BoundingBox(bb.Center - (8 - i) * one18Z, bb.Halfsize - (8 - i) * one18Z);
		bbr[i] = BoundingBox(bb.Center + (i + 1) * one18Z, bb.Halfsize - (i + 1) * one18Z);
	}

	std::vector<Triangle*> leftTriangles[8], rightTriangles[8];

	// add triangles to partitions
	for (Triangle* triangle : triangles)
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


BVHTreeNode* BVHTree::SplitNode(const std::vector<Triangle*>& triangles, unsigned int minTriangles, unsigned int maxDepth, BoundingBox bb)
{
	if (bb.Halfsize.X >= bb.Halfsize.Y)
		if (bb.Halfsize.X > bb.Halfsize.Z)
			return CreateNodeX(triangles, minTriangles, maxDepth);

	if (bb.Halfsize.Y >= bb.Halfsize.Z)
		return CreateNodeY(triangles, minTriangles, maxDepth);
	else
		return CreateNodeZ(triangles, minTriangles, maxDepth);

}

Triangle* BVHTreeNode::Query(const Ray& ray, float& t) const
{
	float tBox = std::numeric_limits<float>::min();
	if (!Intersects(ray, bb, tBox) || tBox > t)
		return nullptr;

	Triangle* triangle = nullptr;

	for (Triangle* tri : triangles)
	{
		float tCurr = 0;
		if (Intersects(ray, tri, tCurr) && tCurr < t)
		{
			t = tCurr;
			triangle = tri;
		}
	}

	if (left != nullptr)
	{
		Triangle* tri = left->Query(ray, t);
		triangle = tri == nullptr ? triangle : tri;
	}

	if (right != nullptr)
	{
		Triangle* tri = right->Query(ray, t);
		triangle = tri == nullptr ? triangle : tri;
	}

	return triangle;
}

Triangle* BVHTree::Query(const Ray& ray, float& t) const
{
	t = std::numeric_limits<float>::max();
	return root->Query(ray, t);
}