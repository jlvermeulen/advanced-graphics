#ifndef BVH_h
#define BVH_h

#include "BBox.h"
#include <vector>
#include <stdint.h>
#include "IntersectionInfo.h"
#include "Ray.h"
#include "Triangle.h"

//! Node descriptor for the flattened tree
struct BVHFlatNode {
  BBox bbox;
  uint32_t start, nPrims, rightOffset;
};

//! \author Brandon Pelfrey
//! A Bounding Volume Hierarchy system for fast Ray-Object intersection tests
class BVH {
  uint32_t nNodes=0, nLeafs=0, leafSize=0;
  std::vector<Triangle*>* build_prims = {};

  //! Build the BVH tree out of build_prims
  void build();

  // Fast Traversal System
  BVHFlatNode *flatTree = nullptr;

  public:
	  BVH(std::vector<Triangle*>* objects, uint32_t leafSize = 4);
  bool getIntersection(const Ray& ray, IntersectionInfo *intersection, bool occlusion) const ;

  ~BVH();
};

#endif
