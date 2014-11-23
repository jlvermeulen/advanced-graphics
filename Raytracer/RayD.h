#pragma once

#include "Vector3D.h"
#include "ColorD.h"
#include "TriangleD.h"

struct RayD
{
public:
	Vector3D Origin, Direction;
	ColorD Color;

	RayD(const Vector3D& origin, const Vector3D& direction, const ColorD& color);

	bool Intersects(const TriangleD& triangle, double& t) const;

private:
};