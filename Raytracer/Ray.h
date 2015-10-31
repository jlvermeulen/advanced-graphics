#pragma once

#include "Vector3F.h"
#include "Color3F.h"
#include "Triangle.h"

struct Ray
{
public:
	Vector3F Origin, Direction, InverseDirection;

	Ray() { }
	Ray(const Vector3F& origin, const Vector3F& direction) : Origin(origin), Direction(direction), InverseDirection(Vector3F(1 / direction.X, 1 / direction.Y, 1 / direction.Z)) {}
	//Ray(const Vector3F& origin, const Vector3F& direction, const Color3F& color) : Origin(origin), Direction(direction), InverseDirection(Vector3F(1 / direction.X, 1 / direction.Y, 1 / direction.Z)), Color(color) { }

	void Reflect(const Vector3F& point, const Vector3F& normal)
	{
		Direction += 2 * Vector3F::Dot(normal, -Direction) * normal;
		Origin = point;
		UpdateInverseDirection();
	}
	Ray Reflect(Ray ray, const Vector3F& point, const Vector3F& normal) { ray.Reflect(point, normal); return ray; }

	void Refract(Vector3F point, const Vector3F& normal, const float& n1, const float& n2)
	{
		float r, c1 = -Vector3F::Dot(Direction, normal), c = c1;
		if (c1 > 0) // entering object
			r = n1 / n2;
		else // leaving object
		{
			r = n2 / n1;
			c = -c1;
		}

		float radicand = 1 - r * r * (1 - c * c);
		if (radicand < 0)
		{
			Reflect(point, c > 0 ? normal : -normal);
			return;
		}

		Direction *= r;
		Direction += (r * c - sqrt(radicand)) * normal;
		Origin = point + (c1 > 0 ? -0.01f : 0.01f) * normal;
		UpdateInverseDirection();
	}

	Ray Refract(Ray ray, const Vector3F& point, const Vector3F& normal, const float& n1, const float& n2) { ray.Refract(point, normal, n1, n2); return ray; }

private:
	void UpdateInverseDirection() { InverseDirection = Vector3F(1 / Direction.X, 1 / Direction.Y, 1 / Direction.Z); }
};