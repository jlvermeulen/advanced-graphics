#pragma once

#include "Material.h"
#include "Octree.h"
#include "Triangle.h"

#include <deque>

struct Object
{
  Object() :
    material()
  { }

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
    octree = Octree(triangles, minTriangles, maxDepth);
  }

  Octree octree;
  Material material;
  std::deque<Triangle> triangles;
};