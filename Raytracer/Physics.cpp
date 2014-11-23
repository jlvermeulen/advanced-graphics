#include "Physics.h"

//---------------------------------------------------------------------------------------
Vector3D Physics::Reflect(const Vector3D& incident, const Vector3D& normal)
{
	return Vector3D(incident - 2 * Vector3D::Dot(incident, normal) * normal);
}

//---------------------------------------------------------------------------------------
Vector3D Physics::Refract(const Vector3D& incident, const Vector3D& normal, const double& n1, const double& n2)
{
	double r = n1 / n2;
	double c = -Vector3D::Dot(normal, incident);
	double c2 = 1 - r * r * (1 - c * c);

	if (c2 < 0)
		return Reflect(incident, normal);

	return Vector3D(r * incident + (r * c - sqrt(c2)) * normal);

}