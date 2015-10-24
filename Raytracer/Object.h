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
	  bvhTree = new BVHTree(triangles, minTriangles, maxDepth);
  }

  BVHTree* bvhTree;
  Octree* octree;
  Material material;
  std::deque<Triangle> triangles;
};