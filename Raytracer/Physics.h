#pragma once

#include <math.h>
#include "Vector3D.h"

class Physics
{
public:
	static Vector3D Reflect(const Vector3D& incident, const Vector3D& normal);
	static Vector3D Refract(const Vector3D& incident, const Vector3D& normal, const double& n1, const double& n2);
};