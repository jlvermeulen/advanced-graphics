#pragma once

#include "Material.h"
#include "Octree.h"
#include "BVHTree.h"
#include "Triangle.h"

#include <deque>

struct Object
{
  Object() :
    material()
  {
  }

  Object(const std::deque<Triangle>& triangles) :
    triangles(triangles),
    material()
  {
  }

  Object(const std::deque<Triangle>& triangles, Material material) :
    triangles(triangles),
    material(material)
  {
  }

  void ConstructOctree(int minTriangles, int maxDepth)
  {
    octree = new Octree(triangles, minTriangles, maxDepth);
  }

  void ConstructBVHtree(int minTriangles, int maxDepth)
  {
	  std::deque<Triangle*> tris;
	  for (int i = 0; i < triangles.size(); i++)
		  tris.push_back(&triangles[i]);
	  bvhTree = new BVHTree(tris, minTriangles, maxDepth);
  }

  BVHTree* bvhTree;
  Octree* octree;
  Material material;
  std::deque<Triangle> triangles;
};