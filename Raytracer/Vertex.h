#pragma once

#include "Vector2F.h"
#include "Vector3F.h"
#include "Color3F.h"

struct Vertex
{
public:
	Vector3F Position, Normal;
	Vector2F UV;

	Vertex() { }
	Vertex(const Vector3F& position, const Vector3F& normal, const Vector2F& uv) : Position(position), Normal(normal), UV(uv) { }
};