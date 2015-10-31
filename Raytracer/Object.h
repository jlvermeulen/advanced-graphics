#pragma once

#include "Material.h"
#include "Octree.h"
#include "Triangle.h"

#include <vector>

struct Object
{
	Object() : material(), octree(nullptr) { }

	Object(const std::vector<Triangle*>& triangles) :
		triangles(triangles),
		material(),
		octree(nullptr)
	{
	}

	Object(const std::vector<Triangle*>& triangles, Material material) :
		triangles(triangles),
		material(material),
		octree(nullptr)
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
			delete octree;

		octree = new Octree(triangles, minTriangles, maxDepth);
	}

	Octree* octree;
	Material material;
	std::vector<Triangle*> triangles;
};