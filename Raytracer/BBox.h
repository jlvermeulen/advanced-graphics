#pragma once

#include "Ray.h"
#include "Vector3.h"
#include <stdint.h>
#include <stdio.h>
__declspec(align(32))struct BBox {
	Vector3 __declspec(align(16))min, max, extent;
  BBox() { }
  BBox(const Vector3& min, const Vector3& max);
  BBox(const Vector3& p);

  bool intersect(const Ray& ray, float *tnear, float *tfar) const;
  void expandToInclude(const Vector3& p);
  void expandToInclude(const BBox& b);
  uint32_t maxDimension() const;
  float surfaceArea() const;
};

