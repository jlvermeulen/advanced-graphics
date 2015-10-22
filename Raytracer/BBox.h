#pragma once

#include "Ray.h"
#include "Vector3.h"
#include <stdint.h>
#include <stdio.h>
__declspec(align(32))struct BBox {
	Vector3 __declspec(align(16))min = Vector3(0.0f, 0.0f, 0.0f), max = Vector3(0.0f, 0.0f, 0.0f), extent = Vector3(0.0f, 0.0f, 0.0f);
  BBox() { }
  BBox(const Vector3& min, const Vector3& max);
  BBox(const Vector3& p);

  bool intersect(const Ray& ray, float *tnear, float *tfar) const;
  void expandToInclude(const Vector3& p);
  void expandToInclude(const BBox& b);
  uint32_t maxDimension() const;
  float surfaceArea() const;
};

