#pragma once

#include "Ray.h"
#include "Triangle.h"
#include "BoundingBox.h"

bool Intersects(const Ray& ray, const Triangle& triangle, double& t);
bool Intersects(const Triangle& triangle, const BoundingBox& boundingBox);