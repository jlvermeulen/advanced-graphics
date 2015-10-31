#pragma once

#include "Ray.h"
#include "Triangle.h"
#include "BoundingBox.h"

#define EPSILON 0.000001

bool Intersects(const Ray& ray, const Triangle& triangle, double& t);
bool Intersects(const Ray& ray, const BoundingBox& bb, double& t);
bool Intersects(const Triangle& triangle, const BoundingBox& boundingBox);