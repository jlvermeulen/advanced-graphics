#pragma once

#include "Vector3D.h"
#include "ColorD.h"

struct RayD
{
public:
	Vector3D Origin, Direction;
	ColorD Color;

	RayD(const Vector3D& origin, const Vector3D& direction, const ColorD& color);

private:
};