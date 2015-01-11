#pragma once

#include "Camera.h"
#include "ColorD.h"
#include "Material.h"
#include "Object.h"
#include "Octree.h"
#include "Triangle.h"

#include <deque>
#include <vector>
#include <utility>

typedef unsigned char uchar;

enum RayDistributionType
{
  none,
  gaussian,
  jitteredStratification,
  stratification,
  uniform
};

struct Intersection
{
  Intersection(const Ray& ray, double time, const Triangle& hit, const Material& material) :
    hitPoint(ray.Origin + time * ray.Direction),
    hit(hit),
    hitMaterial(material)
  {
  }

  Triangle hit;
  Material hitMaterial;
  Vector3D hitPoint;
};

struct Light
{
  Light(const Vector3D& position, const ColorD& color) :
    position(position),
    color(color)
  {
  }

  Vector3D position;
  ColorD color;
};

class Scene
{
public:
  Scene();
  ~Scene();

  bool Render(unsigned char* imageData, bool useOctree, int minTriangles, int maxDepth, int samplesPerPixel, double sigma);
  void LoadDefaultScene();

private:
  ColorD radiance(const Intersection& intersection, Ray ray, double refractiveIndex, int recursionDepth) const;
  ColorD traceRay(Ray ray, double refractiveIndex, int recursionDepth) const;

  void tracePixels(std::pair<ColorD, double>* pixelData, int samplesPerPixel, double sigma);

  double gaussianWeight(double dx, double dy, double sigma) const;

  ColorD calculateDiffuse(const Intersection& intersection) const;
  ColorD calculateReflection(const Intersection& intersection, const Ray& ray, double refractiveIndex, int recursionDepth) const;
  ColorD calculateRefraction(const Intersection& intersection, const Ray& ray, double refractiveIndex, int recursionDepth) const;

public:
  Camera camera;
  std::deque<Object> objects;
  std::vector<Light> lights;
  Object checkerboard;

private:
  bool useOctree_;
};