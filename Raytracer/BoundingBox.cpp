#include "BoundingBox.h"
#include <limits>

BoundingBox::BoundingBox() { }
BoundingBox::BoundingBox(const Vector3D& center, const Vector3D& halfsize) : Center(center), Halfsize(halfsize) { }

BoundingBox BoundingBox::FromTriangles(const std::vector<Triangle*>& triangles)
{
	Vector3D min(std::numeric_limits<double>::max(), std::numeric_limits<double>::max(), std::numeric_limits<double>::max());
	Vector3D max(-std::numeric_limits<double>::max(), -std::numeric_limits<double>::max(), -std::numeric_limits<double>::max());

	for (int i = 0; i < triangles.size(); i++)
		for (int j = 0; j < 3; j++)
		{
			const Vector3D& v = triangles[i]->Vertices[j].Position;
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