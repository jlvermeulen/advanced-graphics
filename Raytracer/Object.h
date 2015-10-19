#pragma once

#include "Material.h"
#include "Octree.h"
#include "Triangle.h"
#include "BVH.h";

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

  void ConstructBVH(int minTriangles)
  {
	  std::vector<Triangle*> bvhTriangles;
	  for (Triangle t : triangles)
		  bvhTriangles.push_back(&t);
	  bvh = new BVH(&bvhTriangles, minTriangles);
  }

  Octree* octree;
  BVH* bvh;
  Material material;
  std::deque<Triangle> triangles;
};