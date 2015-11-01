#pragma once

#include <vector>
#include <limits>
#include "Vector3F.h"
#include "Triangle.h"

struct BoundingSphere
{
public:
	float Radius;
	Vector3F Center;
	BoundingSphere(const Vector3F& center, const float& radius) : Center(center), Radius(radius) { };

	static BoundingSphere FromTriangles(const std::vector<Triangle*>& triangles);
};