#pragma once

#include "Ray.h"
#include "Triangle.h"
#include "BoundingBox.h"

#define EPSILON 0.000002f

bool Intersects(const Ray& ray, const Triangle& triangle, float& t);
bool Intersects(const Ray& ray, const BoundingBox& bb, float& t);
bool Intersects(const Triangle& triangle, const BoundingBox& boundingBox);