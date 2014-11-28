#pragma once

#include <vector>
#include "Vector3D.h"
#include "TriangleD.h"

struct BoundingBox
{
public:
	Vector3D Center, Halfsize;

	BoundingBox() { }
	BoundingBox(const Vector3D& center, const Vector3D& halfsize) : Center(center), Halfsize(halfsize) { }

	static BoundingBox FromTriangles(const std::vector<TriangleD>& triangles)
	{
		Vector3D min(std::numeric_limits<double>::max(), std::numeric_limits<double>::max(), std::numeric_limits<double>::max());
		Vector3D max(std::numeric_limits<double>::min(), std::numeric_limits<double>::min(), std::numeric_limits<double>::min());

		for (std::vector<TriangleD>::const_iterator it = triangles.begin(); it != triangles.end(); ++it)
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

private:
};