#include "Vertex.h"

Vertex::Vertex() :
	Position(),
	Normal(),
	Color(),
	UV()
{ }

Vertex::Vertex(Vector3D position, Vector3D normal, ColorD color, Vector3D uv) :
	Position(position),
	Normal(normal),
	Color(color),
	UV(uv)
{ }

Vertex::Vertex(const Vertex& rhs) :
	Position(rhs.Position),
	Normal(rhs.Normal),
	Color(rhs.Color),
	UV(rhs.UV)
{ }

Vertex::~Vertex()
{
}