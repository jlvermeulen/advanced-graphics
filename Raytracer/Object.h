#pragma once

#include "Material.h"
#include "Octree.h"
#include "Triangle.h"

#include <deque>

struct Object
{
	Object() :
		material()
	{
	}

	Object(const std::deque<Triangle>& triangles) :
		triangles(triangles),
		material()
	{
	}

	Object(const std::deque<Triangle>& triangles, Material* material) :
		triangles(triangles),
		material(material)
	{
	}

	Object(const Object& obj)
	{
		triangles = obj.triangles;
		if (obj.material == nullptr)
			material = nullptr;
		else
			SetMaterial(new Material(*obj.material));
	}

	~Object() { delete material; }

	void SetMaterial(Material* material)
	{
		this->material = material;
		for (Triangle& t : triangles)
			t.Material = material;
	}

	Material* material;
	std::deque<Triangle> triangles;
};