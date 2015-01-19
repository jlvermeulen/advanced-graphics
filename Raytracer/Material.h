#pragma once

#include "ColorD.h"

enum ReflectionType
{
  diffuse,
  specular,
  glossy,
  refractive
};

struct Material
{
  Material() :
    reflType(ReflectionType::diffuse),
    color(1.0, 1.0, 1.0),
    emission(0.0, 0.0, 0.0),
    refrIndex(1.0),
    specularExponent(0.0),
    transparency(0.0)
  {
  }

  Material(ReflectionType reflType, const ColorD& color, const ColorD& emission, double refrIndex, double specularExponent, double transparency) :
    reflType(reflType),
    color(color),
    emission(emission),
    refrIndex(refrIndex),
    specularExponent(specularExponent),
    transparency(transparency)
  {
  }

  ReflectionType reflType;
  ColorD color;
  ColorD emission;
  double refrIndex;
  double specularExponent;
  double transparency;
};