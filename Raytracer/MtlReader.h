#pragma once

#include "ColorD.h"
#include "Material.h"
#include "Reader.h"

#include <map>
#include <vector>

namespace MtlType
{
  enum MtlType
  {
    ambient,  // Emission
    diffuse,  // Color
    dissolve,
    illumination,
    name,
    refractiveIndex,
    transparency,
    nil
  };
}

class MtlReader : public Reader
{
public:
  MtlReader();
  ~MtlReader();

  std::map<std::string, Material> parseFile(const char* fileName);
  void parseLine(const std::vector<std::string>& segments);

private:
  void parseAmbient(CVSIterator& it);
  void parseDiffuse(CVSIterator& it);
  void parseDissolve(CVSIterator& it);
  void parseIllumination(CVSIterator& it);
  void parseName(CVSIterator& it);
  void parseRefractiveIndex(CVSIterator& it);
  void parseTransparency(CVSIterator& it);

  ColorD parseColor(CVSIterator& it) const;
  MtlType::MtlType parseType(const std::string& type) const;

private:
  std::string current_;
  std::map<std::string, Material> materials_;
};