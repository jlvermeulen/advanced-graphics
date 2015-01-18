#pragma once

#include "Camera.h"
#include "ColorD.h"
#include <deque>
#include "Octree.h"
#include "Triangle.h"
#include <vector>
#include <utility>
#include <tuple>
#include <random>

typedef unsigned char uchar;

enum ReflectionType
{
  diffuse,
  specular,
  refractive
};

enum RayDistributionType
{
  none,
  gaussian,
  jitteredStratification,
  stratification,
  uniform
};

struct Material
{
  Material() :
    reflType(ReflectionType::diffuse),
    color(1.0, 1.0, 1.0),
    emission(0.0, 0.0, 0.0),
    refrIndex(1.0),
    transparency(0.0)
  {
  }

  Material(const ReflectionType& reflType, const ColorD& color, const ColorD& emission, float refrIndex, float transparency) :
    reflType(reflType),
    color(color),
    emission(emission),
    refrIndex(refrIndex),
    transparency(transparency)
  {
  }

  ReflectionType reflType;
  ColorD color;
  ColorD emission;
  float refrIndex;
  float transparency;
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
struct Lightarea
{
	Lightarea(const Triangle& triangle, const ColorD& color, const float& emission) :
		triangle(triangle),
		color(color),
		emission(emission)
	{
	}

	Triangle triangle;
	ColorD color;
	float emission;
};
struct Object
{
  Object() : material(ReflectionType::diffuse, ColorD(), ColorD(), 1, 0) { }

  Object(const std::deque<Triangle>& triangles, Material material) :
    triangles(triangles),
    material(material)
  {
  }

  void CreateOctree(int minTriangles, int maxDepth)
  {
    octree = new Octree(triangles, minTriangles, maxDepth);
  }

  Octree* octree;
  Material material;
  std::deque<Triangle> triangles;

};

class Scene
{
public:
  Scene();
  ~Scene();

  bool Render(unsigned char* imageData, int minTriangles, int maxDepth, int samplesPerPixel, double sigma);
  void LoadDefaultScene();

private:
  ColorD TraceRay(const Ray& ray);
  ColorD ComputeRadiance(const Vector3D& point, const Vector3D& in, const Triangle& triangle, const Material& material, unsigned int depth);
  ColorD DirectIllumination(const Vector3D& point, const Material& material);
  ColorD IndirectIllumination(const Vector3D& point, const Vector3D& in, const Triangle& triangle, const Material& material, unsigned int depth);

  void TracePixels(std::pair<ColorD, double>* pixelData, int samplesPerPixel, double sigma);

  double GaussianWeight(double dx, double dy, double sigma) const;

  /*ColorD calculateDiffuse(const Intersection& intersection) const;
  ColorD calculateReflection(const Intersection& intersection, const Ray& ray, double refractiveIndex, int recursionDepth) const;
  ColorD calculateRefraction(const Intersection& intersection, const Ray& ray, double refractiveIndex, int recursionDepth) const;*/
  std::tuple<Ray, Triangle, double> SampleLight(const Vector3D& hitPoint);
  bool Scene::FirstHitInfo(const Ray& ray, double& time, Triangle& triangle, Material& mat) const;

public:
  Camera camera;
  std::vector<Object> objects;
  //std::vector<Light> lights;
  Object checkerboard;
  //std::vector<Lightarea> lightareas;
  std::deque<Object> lights;
  std::uniform_real_distribution<double> dist;
  std::mt19937 gen;

private:
};