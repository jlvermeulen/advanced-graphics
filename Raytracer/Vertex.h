#pragma once

#include "Vector2D.h"
#include "Vector3D.h"
#include "ColorD.h"

struct Vertex
{
public:
	Vector3D Position, Normal;
	Vector2D UV;

	Vertex() { }
	Vertex(const Vector3D& position, const Vector3D& normal, const Vector2D& uv) : Position(position), Normal(normal), UV(uv) { }
};