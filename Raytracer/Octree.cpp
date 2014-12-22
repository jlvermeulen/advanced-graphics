#include "Octree.h"
#include "Intersections.h"

OctreeNode::OctreeNode() { }
OctreeNode::OctreeNode(const std::vector<Triangle>& triangles, const BoundingBox& bb, unsigned int minTriangles, unsigned int maxDepth)
{
	this->bb = bb;

	if (triangles.size() <= minTriangles || maxDepth == 0)
	{
		this->triangles = triangles;
		this->children = nullptr;
		return;
	}

	children = new OctreeNode[8];
	Vector3D half = bb.Halfsize / 2;
	for (int i = 0; i < 8; i++)
	{
		double x = bb.Center.X + (i & 4 ? half.X : -half.X);
		double y = bb.Center.Y + (i & 2 ? half.Y : -half.Y);
		double z = bb.Center.Z + (i & 1 ? half.Z : -half.Z);
		BoundingBox childBB(Vector3D(x, y, z), half);

		std::vector<Triangle> childTriangles;
		for(std::vector<Triangle>::const_iterator it = triangles.begin(); it != triangles.end(); ++it)
			if (Intersects(*it, childBB))
				childTriangles.push_back(*it);

		children[i] = OctreeNode(childTriangles, childBB, minTriangles, maxDepth - 1);
	}
}

OctreeNode::~OctreeNode() { delete [] children; }

bool OctreeNode::Query(const Ray& ray, Triangle& triangle, double& t) const
{
	double tMin;
	if (!Intersects(ray, bb, tMin))
		return false;

	bool hit = false;
	const Triangle* tri = nullptr;
	for(std::vector<Triangle>::const_iterator it = triangles.begin(); it != triangles.end(); ++it)
	{
		double tCurr = 0;
		if (Intersects(ray, *it, tCurr) && tCurr < tMin)
		{
			tMin = tCurr;
			tri = &*it;
			hit = true;
		}
	}

	if (children != nullptr)
		for (int i = 0; i < 8; i++)
			hit &= children[i].Query(ray, triangle, tMin);

	triangle = *tri;
	t = tMin;
	return hit;
}

Octree::Octree() { }
Octree::Octree(const std::vector<Triangle>& triangles, int minTriangles, int maxDepth) : root(OctreeNode(triangles, BoundingBox::FromTriangles(triangles), minTriangles, maxDepth)) { }

bool Octree::Query(const Ray& ray, Triangle& triangle, double& t) const { return root.Query(ray, triangle, t); }