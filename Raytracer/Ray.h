#pragma once

#include "Vector3D.h"
#include "ColorD.h"
#include "Triangle.h"

struct Ray
{
public:
	Vector3D Origin, Direction;
	ColorD Color;

	Ray(const Vector3D& origin, const Vector3D& direction, const ColorD& color) : Origin(origin), Direction(direction), Color(color) { }

private:
};