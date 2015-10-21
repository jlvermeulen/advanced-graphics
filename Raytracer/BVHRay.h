#pragma once

#include "Vector3.h"

__declspec(align(16))struct BVHRay {
	__declspec(align(16))Vector3 o; // Ray Origin
	__declspec(align(16))Vector3 d; // Ray Direction
	__declspec(align(16))Vector3 inv_d; // Inverse of each Ray Direction component

	BVHRay(const Vector3& o, const Vector3& d)
		: o(o), d(d), inv_d(Vector3(1, 1, 1).cdiv(d)) { }
};
