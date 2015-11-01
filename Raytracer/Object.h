#pragma once

#include "Material.h"
#include "Octree.h"
#include "BVHTree.h"
#include "Triangle.h"

#include <vector>

struct Object
{
	Object() : material(), octree(nullptr), bvhTree(nullptr) { }

	Object(const std::vector<Triangle*>& triangles) :
		triangles(triangles),
		material(),
		octree(nullptr),
		bvhTree(nullptr)
	{
	}

	Object(const std::vector<Triangle*>& triangles, Material material) :
		triangles(triangles),
		material(material),
		octree(nullptr),
		bvhTree(nullptr)
	{
	}

	~Object()
	{
		for (unsigned int i = 0; i < triangles.size(); i++)
			delete triangles[i];

		if (octree != nullptr)
			delete octree;
	}

	void ConstructOctree(int minTriangles, int maxDepth)
	{
		if (octree != nullptr)
			return;

		octree = new Octree(triangles, minTriangles, maxDepth);
	}

	void ConstructBVHTree(int minTriangles, int maxDepth)
	{
		if (bvhTree != nullptr)
			return;

		bvhTree = new BVHTree(triangles, minTriangles, maxDepth);
	}

	BVHTree* bvhTree;
	Octree* octree;
	Material material;
	std::vector<Triangle*> triangles;
};