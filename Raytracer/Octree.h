#pragma once

#include "TriangleD.h"
#include "BoundingBox.h"

struct OctreeNode
{
public:
	TriangleD* Triangles;
	int NTriangles;
	OctreeNode* Children;
	BoundingBox BoundingBox;

	OctreeNode();
	OctreeNode(TriangleD* triangles, int nTriangles);

private:
};

struct Octree
{
public:
	Octree();
	Octree(TriangleD* triangles, int nTriangles);

private:
	OctreeNode root;
};