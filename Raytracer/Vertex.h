#pragma once

#include "Vector3D.h"
#include "ColorD.h"

struct Vertex
{
public:
	Vector3D Position, Normal;
	Vector3D UV; // TODO: make Vector2D

	Vertex() { }
	Vertex(const Vector3D& position, const Vector3D& normal, const ColorD& color, const Vector3D& uv) : Position(position), Normal(normal), UV(uv) { }

private:
};