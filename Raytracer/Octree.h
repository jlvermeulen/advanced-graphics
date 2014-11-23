#pragma once

#include "TriangleD.h"

struct OctreeNode
{
public:
	TriangleD* Triangles;
	int NTriangles;
	OctreeNode* Children;

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