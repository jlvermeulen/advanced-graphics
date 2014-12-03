#pragma once

#include <vector>
#include "Vector3D.h"
#include "Triangle.h"

struct BoundingBox
{
public:
	Vector3D Center, Halfsize;

	BoundingBox();
	BoundingBox(const Vector3D& center, const Vector3D& halfsize) : Center(center), Halfsize(halfsize);

	static BoundingBox FromTriangles(const std::vector<Triangle>& triangles);

private:
};