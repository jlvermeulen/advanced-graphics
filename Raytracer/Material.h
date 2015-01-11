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
    refrIndex(1.0),
    transparency(0.0)
  {
  }

  Material(ReflectionType reflType, const ColorD& ambient, const ColorD& diffuse, const ColorD& specular, double specularWeight, double refrIndex, double transparency) :
    reflType(reflType),
    ambient(ambient),
    diffuse(diffuse),
    specular(specular),
    specularWeight(specularWeight),
    refrIndex(refrIndex),
    transparency(transparency)
  {
  }

  ReflectionType reflType;
  ColorD ambient;
  ColorD diffuse;
  ColorD specular;
  double specularWeight;
  double refrIndex;
  double transparency;
};