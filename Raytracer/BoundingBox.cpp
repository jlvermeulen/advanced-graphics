#include "BoundingBox.h"
#include <limits>

BoundingBox::BoundingBox() { }
BoundingBox::BoundingBox(const Vector3F& center, const Vector3F& halfsize) : Center(center), Halfsize(halfsize) { }

BoundingBox BoundingBox::FromTriangles(const std::vector<Triangle*>& triangles)
{
	Vector3F min(std::numeric_limits<float>::max(), std::numeric_limits<float>::max(), std::numeric_limits<float>::max());
	Vector3F max(-std::numeric_limits<float>::max(), -std::numeric_limits<float>::max(), -std::numeric_limits<float>::max());

	for (unsigned int i = 0; i < triangles.size(); i++)
		for (unsigned int j = 0; j < 3; j++)
		{
			const Vector3F& v = triangles[i]->Vertices[j].Position;
			for (int j = 0; j < 3; j++)
			{
				if (v[j] < min[j])
					min[j] = v[j];
				if (v[j] > max[j])
					max[j] = v[j];
			}
		}

	return BoundingBox((max + min) / 2, (max - min) / 2);
}