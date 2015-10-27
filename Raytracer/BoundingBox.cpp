#include "BoundingBox.h"
#include <limits>

BoundingBox::BoundingBox() { }
BoundingBox::BoundingBox(const Vector3D& center, const Vector3D& halfsize) : Center(center), Halfsize(halfsize) { }

BoundingBox BoundingBox::FromTriangles(const std::deque<Triangle>& triangles)
{
	Vector3D min(std::numeric_limits<double>::max(), std::numeric_limits<double>::max(), std::numeric_limits<double>::max());
	Vector3D max(-std::numeric_limits<double>::max(), -std::numeric_limits<double>::max(), -std::numeric_limits<double>::max());


	for (std::deque<Triangle>::const_iterator it = triangles.begin(); it != triangles.end(); ++it)
		for (int i = 0; i < 3; i++)
		{
			const Vector3D& v = it->Vertices[i].Position;
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

double BoundingBox::SurfaceArea()
{
	double xedge = 2 * Halfsize.X;
	double yedge = 2 * Halfsize.Y;
	double zedge = 2 * Halfsize.Z;

	return 2 * (xedge*yedge + xedge*zedge + yedge*zedge);
}