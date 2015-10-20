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
	  bvhTriangles.clear();
	  for (int i = 0; i < triangles.size();++i)
		  bvhTriangles.push_back(&triangles[i]);
	  bvh = new BVH(&bvhTriangles, minTriangles);
  }
  std::vector<Triangle*> bvhTriangles;
  Octree* octree;
  BVH* bvh;
  Material material;
  std::deque<Triangle> triangles;
};