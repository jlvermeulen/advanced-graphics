#pragma once

#include "RayD.h"
#include "TriangleD.h"
#include "BoundingBox.h"

bool Intersects(const RayD& ray, const TriangleD& triangle, double& t);
bool Intersects(const TriangleD& triangle, const BoundingBox& boundingBox);