#pragma once

#include <vector>
#include "TriangleD.h"
#include "BoundingBox.h"

class OctreeNode
{
public:
	OctreeNode();
	OctreeNode(std::vector<TriangleD> triangles);

private:
	std::vector<TriangleD> Triangles;
	OctreeNode* Children;
	BoundingBox BoundingBox;
};

class Octree
{
public:
	Octree();
	Octree(std::vector<TriangleD> triangles);

private:
	OctreeNode root;
};