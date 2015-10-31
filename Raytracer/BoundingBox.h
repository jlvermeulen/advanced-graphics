#pragma once

#include <vector>
#include <limits>
#include "Vector3F.h"
#include "Triangle.h"

struct BoundingBox
{
public:
	Vector3F Center, Halfsize;

	BoundingBox();
	BoundingBox(const Vector3F& center, const Vector3F& halfsize);

	static BoundingBox FromTriangles(const std::vector<Triangle*> & triangles);

private:
};
