#pragma once

#include "Vector3D.h"
#include "ColorD.h"
#include "Triangle.h"

struct Ray
{
public:
	Vector3D Origin, Direction, InverseDirection;

	Ray() { }
	Ray(const Vector3D& origin, const Vector3D& direction) : Origin(origin), Direction(direction), InverseDirection(Vector3D(1 / direction.X, 1 / direction.Y, 1 / direction.Z)) {}
	//Ray(const Vector3D& origin, const Vector3D& direction, const ColorD& color) : Origin(origin), Direction(direction), InverseDirection(Vector3D(1 / direction.X, 1 / direction.Y, 1 / direction.Z)), Color(color) { }

	void Reflect(const Vector3D& point, const Vector3D& normal)
	{
		Direction += 2 * Vector3D::Dot(normal, -Direction) * normal;
		Origin = point;
		UpdateInverseDirection();
	}
	Ray Reflect(Ray ray, const Vector3D& point, const Vector3D& normal) { ray.Reflect(point, normal); return ray; }

	void Refract(Vector3D point, const Vector3D& normal, const double& n1, const double& n2)
	{
		double r, c1 = -Vector3D::Dot(Direction, normal), c = c1;
		if (c1 > 0) // entering object
			r = n1 / n2;
		else // leaving object
		{
			r = n2 / n1;
			c = -c1;
		}

		double radicand = 1 - r * r * (1 - c * c);
		if (radicand < 0)
		{
			Reflect(point, c > 0 ? normal : -normal);
			return;
		}

		Direction *= r;
		Direction += (r * c - sqrt(radicand)) * normal;
		Origin = point + (c1 > 0 ? -0.01 : 0.01) * normal;
		UpdateInverseDirection();
	}

	Ray Refract(Ray ray, const Vector3D& point, const Vector3D& normal, const double& n1, const double& n2) { ray.Refract(point, normal, n1, n2); return ray; }

private:
	void UpdateInverseDirection() { InverseDirection = Vector3D(1 / Direction.X, 1 / Direction.Y, 1 / Direction.Z); }
};