#ifndef IntersectionInfo_h_
#define IntersectionInfo_h_

class Triangle;

struct IntersectionInfo {
  float t; // Intersection distance along the ray
  Triangle* triangle; // Object that was hit
  Vector3D hit; // Location of the intersection
};

#endif
