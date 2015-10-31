#pragma once

#include "Color3F.h"

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
    color(1.0f, 1.0f, 1.0f),
    emission(0.0f, 0.0f, 0.0f),
    refrIndex(1.0f),
    specularExponent(0.0f),
    transparency(0.0f)
  {
  }

  Material(ReflectionType reflType, const Color3F& color, const Color3F& emission, float refrIndex, float specularExponent, float transparency) :
    reflType(reflType),
    color(color),
    emission(emission),
    refrIndex(refrIndex),
    specularExponent(specularExponent),
    transparency(transparency)
  {
  }

  ReflectionType reflType;
  Color3F color;
  Color3F emission;
  float refrIndex;
  float specularExponent;
  float transparency;
};