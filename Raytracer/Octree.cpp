#include "Octree.h"
#include "Intersections.h"

OctreeNode::OctreeNode() { }
OctreeNode::OctreeNode(const std::deque<Triangle>& triangles, const BoundingBox& bb, unsigned int minTriangles, unsigned int maxDepth)
{
	this->bb = bb;
	this->children = std::vector<OctreeNode>();

	if (triangles.size() <= minTriangles || maxDepth == 0)
	{
		this->triangles = triangles;
		return;
	}

	Vector3D half = bb.Halfsize / 2;
	for (int i = 0; i < 8; i++)
	{
		double x = bb.Center.X + (i & 4 ? half.X : -half.X);
		double y = bb.Center.Y + (i & 2 ? half.Y : -half.Y);
		double z = bb.Center.Z + (i & 1 ? half.Z : -half.Z);
		BoundingBox childBB(Vector3D(x, y, z), half);

		std::deque<Triangle> childTriangles;
		for(std::deque<Triangle>::const_iterator it = triangles.begin(); it != triangles.end(); ++it)
			if (Intersects(*it, childBB))
				childTriangles.push_back(*it);

		children.push_back(OctreeNode(childTriangles, childBB, minTriangles, maxDepth - 1));
	}
}

bool OctreeNode::Query(const Ray& ray, Triangle& triangle, double& t) const
{
	double tMin;
	if (!Intersects(ray, bb, tMin))
		return false;

	bool hit = false;
	const Triangle* tri = nullptr;
	for(std::deque<Triangle>::const_iterator it = triangles.begin(); it != triangles.end(); ++it)
	{
		double tCurr = 0;
		if (Intersects(ray, *it, tCurr) && tCurr < tMin)
		{
			tMin = tCurr;
			tri = &*it;
			hit = true;
		}
	}

	if (!children.empty())
		for (int i = 0; i < 8; i++)
			hit &= children[i].Query(ray, triangle, tMin);

	triangle = *tri;
	t = tMin;
	return hit;
}

Octree::Octree() { }
Octree::Octree(const std::deque<Triangle>& triangles, int minTriangles, int maxDepth) : root(triangles, BoundingBox::FromTriangles(triangles), minTriangles, maxDepth) { }

bool Octree::Query(const Ray& ray, Triangle& triangle, double& t) const { return root.Query(ray, triangle, t); }