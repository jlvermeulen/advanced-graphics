#pragma once

#include "Vector3D.h"
#include "ColorD.h"

struct Vertex
{
public:
	Vector3D Position, Normal;
	ColorD Color;
	Vector3D UV; // TODO: make Vector2D

private:
};