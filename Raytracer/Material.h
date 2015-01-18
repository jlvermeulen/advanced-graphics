#pragma once

#include "ColorD.h"

enum ReflectionType
{
  diffuse,
  specular,
  refractive
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

  Material(ReflectionType reflType, const ColorD& color, const ColorD& emission, double refrIndex, double transparency) :
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
  double refrIndex;
  double transparency;
};