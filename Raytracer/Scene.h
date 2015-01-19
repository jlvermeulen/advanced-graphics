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
#include <tuple>
#include <random>

typedef unsigned char uchar;

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

  bool Render(unsigned char* imageData, int minTriangles, int maxDepth, int samplesPerPixel, double sigma);
  void LoadDefaultScene();
  void LoadDefaultScene2();

private:
  ColorD TraceRay(const Ray& ray);
  ColorD ComputeRadiance(const Vector3D& point, const Vector3D& in, const Triangle& triangle, const Material& material, unsigned int depth);
  ColorD DirectIllumination(const Vector3D& point, const Vector3D& normal, const Material& material);
  ColorD IndirectIllumination(Vector3D point, const Vector3D& in, const Vector3D& triangle, const Material& material, unsigned int depth);

  void TracePixels(std::pair<ColorD, double>* pixelData, int samplesPerPixel, double sigma);

  double GaussianWeight(double dx, double dy, double sigma) const;

  /*ColorD calculateDiffuse(const Intersection& intersection) const;
  ColorD calculateReflection(const Intersection& intersection, const Ray& ray, double refractiveIndex, int recursionDepth) const;
  ColorD calculateRefraction(const Intersection& intersection, const Ray& ray, double refractiveIndex, int recursionDepth) const;*/
  std::pair<Ray, Triangle> SampleLight(const Vector3D& hitPoint);
  bool Scene::FirstHitInfo(const Ray& ray, double& time, Triangle& triangle, Material& mat) const;

public:
  Camera camera;
  std::deque<Object> objects;
  Object checkerboard;
  std::deque<Object> lights;
  std::uniform_real_distribution<double> dist;
  std::mt19937 gen;

private:
};