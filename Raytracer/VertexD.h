#pragma once

#include "Vector3D.h"
#include "ColorD.h"

struct VertexD
{
public:
	Vector3D Position, Normal;
	ColorD Color;
	Vector3D UV; // TODO: make Vector2D

	VertexD() { }
	VertexD(const Vector3D& position, const Vector3D& normal, const ColorD& color, const Vector3D& uv) : Position(position), Normal(normal), Color(color), UV(uv) { }

private:
};