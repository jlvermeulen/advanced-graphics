#pragma once
#include <deque>
#include <limits>
#include "Vector3D.h"
#include "Triangle.h"
struct BoundingSphere
{
public:
	double radius;
	Vector3D center;
	BoundingSphere();
	~BoundingSphere();
	BoundingSphere(const Vector3D& center, const double& radius) :center(center), radius(radius){};
	static BoundingSphere FromTriangles(const std::deque<Triangle>& triangles);
};

