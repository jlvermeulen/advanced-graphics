#pragma once

#include "Vector3D.h"
#include "ColorD.h"

struct Vertex
{
public:
	Vector3D Position, Normal;
	ColorD Color;
	Vector3D UV; // TODO: make Vector2D

	Vertex();
	Vertex(Vector3D position, Vector3D normal, ColorD color, Vector3D uv);
	Vertex(const Vertex& rhs);
	~Vertex();

private:
};