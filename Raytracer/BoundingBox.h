#pragma once

#include "Vector3D.h"

struct BoundingBox
{
public:
	Vector3D Center, Halfsize;

	BoundingBox() { }
	BoundingBox(const Vector3D& center, const Vector3D& halfsize) : Center(center), Halfsize(halfsize) { }

private:
};