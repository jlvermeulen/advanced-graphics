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

		if (bvhTree != nullptr)
			delete bvhTree;
	}

	void ConstructOctree()
	{
		if (octree != nullptr)
			return;

		octree = new Octree(triangles);
	}

	void ConstructBVHTree()
	{
		if (bvhTree != nullptr)
			return;

		bvhTree = new BVHTree(triangles);
	}

	BVHTree* bvhTree;
	Octree* octree;
	Material material;
	std::vector<Triangle*> triangles;
};