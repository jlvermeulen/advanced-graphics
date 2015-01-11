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

    case MtlType::specular:
      parseSpecular(++it);
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
  materials_.at(current_).ambient = parseColor(it);
}

//--------------------------------------------------------------------------------
void MtlReader::parseDiffuse(CVSIterator& it)
{
  materials_.at(current_).diffuse = parseColor(it);
}

//--------------------------------------------------------------------------------
void MtlReader::parseDissolve(CVSIterator& it)
{
  materials_.at(current_).transparency = 1.0 - parseDouble(it);
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
void MtlReader::parseSpecular(CVSIterator& it)
{
  materials_.at(current_).specular = parseColor(it);
}

//--------------------------------------------------------------------------------
void MtlReader::parseSpecularWeight(CVSIterator& it)
{
  materials_.at(current_).specularWeight = parseDouble(it);
}

//--------------------------------------------------------------------------------
void MtlReader::parseTransparency(CVSIterator& it)
{
  materials_.at(current_).transparency = parseDouble(it);
}

//--------------------------------------------------------------------------------
ColorD MtlReader::parseColor(CVSIterator& it) const
{
  ColorD color;
  
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
    return MtlType::ambient;
  else if (type == "Kd")
    return MtlType::diffuse;
  else if (type == "Ks")
    return MtlType::specular;
  else if (type == "Ns")
    return MtlType::specularWeight;
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