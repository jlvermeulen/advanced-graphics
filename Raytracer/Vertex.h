#pragma once

#include "Vector3D.h"
#include "Color.h"

struct Vertex
{
public:
	Vector3D Position, Normal;
	Color Color;
	Vector3D UV; // TODO: make Vector2D

private:
};