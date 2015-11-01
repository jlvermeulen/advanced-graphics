#include "MtlReader.h"

#include <algorithm>
#include <fstream>
#include <string>

//--------------------------------------------------------------------------------
MtlReader::MtlReader() :
  Reader()
{

}

//--------------------------------------------------------------------------------
MtlReader::~MtlReader()
{

}

//--------------------------------------------------------------------------------
std::map<std::string, Material> MtlReader::parseFile(const char* fileName)
{
  // Open file stream
  std::ifstream fin;
  fin.open(fileName);

  std::string line;

  while (std::getline(fin, line))
  {
    std::vector<std::string> segments = split(line, ' ');

    parseLine(segments);
  }

  return materials_;
}

//--------------------------------------------------------------------------------
void MtlReader::parseLine(const std::vector<std::string>& segments)
{
  CVSIterator it = segments.cbegin();

  // No segments at all
  if (it == segments.cend())
    return;

  switch (parseType(*it))
  {
    case MtlType::ambient:
      parseAmbient(++it);
      break;

    case MtlType::diffuse:
      parseDiffuse(++it);
      break;

    case MtlType::dissolve:
      parseDissolve(++it);
      break;

    case MtlType::illumination:
      parseIllumination(++it);
      break;

    case MtlType::name:
      parseName(++it);
      break;

    case MtlType::refractiveIndex:
      parseRefractiveIndex(++it);
      break;

    case MtlType::specularWeight:
      parseSpecularWeight(++it);
      break;

    case MtlType::transparency:
      parseTransparency(++it);
      break;

    case MtlType::nil:
    default:
      // Don't do anything
      break;
  }
}

//--------------------------------------------------------------------------------
void MtlReader::parseAmbient(CVSIterator& it)
{
  materials_.at(current_).emission = parseColor(it);
}

//--------------------------------------------------------------------------------
void MtlReader::parseDiffuse(CVSIterator& it)
{
  materials_.at(current_).color = parseColor(it);
}

//--------------------------------------------------------------------------------
void MtlReader::parseDissolve(CVSIterator& it)
{
  materials_.at(current_).transparency = 1.0f - parseDouble(it);
}

//--------------------------------------------------------------------------------
void MtlReader::parseIllumination(CVSIterator& it)
{
  int illum = parseInteger(it);
  ReflectionType reflType;

  switch (illum)
  {
    case 1:
    case 2:
    {
      reflType = ReflectionType::diffuse;
      break;
    }

    case 3:
    {
      reflType = ReflectionType::specular;
      break;
    }

    default:
    {
      reflType = ReflectionType::refractive;
      break;
    }
  }

  materials_.at(current_).reflType = reflType;
}

//--------------------------------------------------------------------------------
void MtlReader::parseName(CVSIterator& it)
{
  std::string name = *it;

  materials_.emplace(name, Material());
  current_ = name;
}

//--------------------------------------------------------------------------------
void MtlReader::parseRefractiveIndex(CVSIterator& it)
{
  materials_.at(current_).refrIndex = parseDouble(it);
}

//--------------------------------------------------------------------------------
void MtlReader::parseSpecularWeight(CVSIterator& it)
{
  materials_.at(current_).specularExponent = parseDouble(it);
}

//--------------------------------------------------------------------------------
void MtlReader::parseTransparency(CVSIterator& it)
{
  materials_.at(current_).transparency = parseDouble(it);
}

//--------------------------------------------------------------------------------
Color3F MtlReader::parseColor(CVSIterator& it) const
{
  Color3F color;
  
  color.R = parseDouble(it);
  color.G = parseDouble(++it);
  color.B = parseDouble(++it);

  return color;
}

//--------------------------------------------------------------------------------
MtlType::MtlType MtlReader::parseType(const std::string& type) const
{
  if (type == "newmtl")
    return MtlType::name;
  else if (type == "Ka")
    return MtlType::ambient;  // Emission
  else if (type == "Kd")
    return MtlType::diffuse;  // Color
  else if (type == "Ns")
    return MtlType::specularWeight; // Specular Exponent
  else if (type == "Ni")
    return MtlType::refractiveIndex;
  else if (type == "d")
    return MtlType::dissolve;
  else if (type == "Tr")
    return MtlType::transparency;
  else if (type == "illum")
    return MtlType::illumination;
  else
    return MtlType::nil;
}