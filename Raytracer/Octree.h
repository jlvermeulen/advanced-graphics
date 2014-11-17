#pragma once

#include <vector>
#include "Triangle.h"

using namespace std;

class Octree
{
public:
	Octree(vector<Triangle> triangles);
	~Octree();
};