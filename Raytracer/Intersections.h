#pragma once

#include "Ray.h"
#include "Triangle.h"
#include "BoundingBox.h"
#include "BoundingSphere.h"

bool Intersects(const Ray& ray, const Triangle& triangle, double& t);
bool Intersects(const Ray& ray, const BoundingBox& bb, double& t);
bool Intersects(const Ray& ray, const BoundingSphere& bb, double& t);
bool Intersects(const Triangle& triangle, const BoundingBox& boundingBox);
bool Intersects(const Triangle& triangle, const BoundingSphere& boundingBox);