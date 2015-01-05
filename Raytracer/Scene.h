#pragma once

#include "Camera.h"
#include "ColorD.h"
#include <deque>
#include "Octree.h"
#include "Triangle.h"
#include <vector>

typedef unsigned char uchar;


enum ReflectionType
{
  diffuse,
  specular,
  refractive
};

struct Intersection
{
  Intersection(const Ray& ray, double time, const Triangle& hit) :
    hitPoint(ray.Origin + time * ray.Direction),
    hit(hit)
  {
  }

  Vector3D hitPoint;
  Triangle hit;
};

struct Material
{
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

struct Object
{
  Object(const std::deque<Triangle>& triangles, Material material) :
    triangles(triangles),
    material(material)
  {
  }

  void CreateOctree(int minTriangles, int maxDepth)
  {
    octree = Octree(triangles, minTriangles, maxDepth);
  }

  Octree octree;
  Material material;
  std::deque<Triangle> triangles;

};

class Scene
{
public:
  Scene();
  ~Scene();

  bool Render(unsigned char* imageData, bool useOctree, int minTriangles, int maxDepth);
  void LoadDefaultScene();

private:
  ColorD radiance(const Intersection& intersection, Ray ray, double refractiveIndex, int recursionDepth) const;
  ColorD traceRay(Ray ray, double refractiveIndex, int recursionDepth) const;

  ColorD calculateDiffuse(const Intersection& intersection) const;
  ColorD calculateReflection() const;
  ColorD calculateRefraction() const;

public:
  Camera camera;
  std::vector<Object> objects;
  std::vector<Light> lights;

private:
  bool useOctree_;
};