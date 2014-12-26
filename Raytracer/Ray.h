#pragma once

#include "Vector3D.h"
#include "ColorD.h"
#include "Triangle.h"

struct Ray
{
public:
	Vector3D Origin, Direction, InverseDirection;
	ColorD Color;

	Ray() { }
	Ray(const Vector3D& origin, const Vector3D& direction, const ColorD& color) : Origin(origin), Direction(direction), InverseDirection(Vector3D(1 / direction.X, 1 / direction.Y, 1 / direction.Z)), Color(color) { }

private:
};