#pragma once

#include "Camera.h"
#include "ColorD.h"
#include <deque>
#include "Triangle.h"
#include <vector>

enum ReflectionType
{
  diffuse,
  specular,
  refractive
};

struct Material
{
  Material(const ReflectionType& reflType, const ColorD& color, const ColorD& emission, float refrIndex, float transparency) :
    ReflType(reflType),
    Color(color),
    Emission(emission),
    RefrIndex(refrIndex),
    Transparency(transparency)
  {
  }

  ReflectionType ReflType;
  ColorD Color;
  ColorD Emission;
  float RefrIndex;
  float Transparency;
};

struct Light
{
  Light(const Vector3D& position, const ColorD& color) :
    Position(position),
    Color(color)
  {
  }

  Vector3D Position;
  ColorD Color;
};

struct Object
{
  Object(const std::deque<Triangle>& triangles, Material material) :
    Triangles(triangles),
    Material(material)
  {
  }

  Material Material;
  std::deque<Triangle> Triangles;
};

class Scene
{
public:
  Scene();

public:
  Camera Camera;
  std::vector<Object> Objects;
  std::vector<Light> Lights;
};