#pragma once

#include "Reader.h"
#include "Object.h"
#include "Triangle.h"
#include "Vector3D.h"

#include <map>
#include <vector>

typedef std::vector<Vector3D> Vector3List;

namespace ObjType
{
  enum ObjType
  {
    face,
    normal,
    space,
    texCoords,
    vertex,
    library,
    material,
    nil
  };
}

class ObjReader : public Reader
{
public:
  ObjReader();
  ~ObjReader();

public:
  std::vector<Object*> parseFile(const char* fileName);
  void parseLine(const std::string& path, const std::vector<std::string>& segments);

  void parseLibrary(const std::string& path, CVSIterator& it);
  void parseMaterial(CVSIterator& it);
  void parseFace(CVSIterator& it);
  void parseNormal(CVSIterator& it);
  void parseTexCoords(CVSIterator& it);
  void parseVertex(CVSIterator& it);

  ObjType::ObjType parseType(const std::string& type) const;

private:
  void clearObject();
  void reset();

private:
  std::vector<Vector3D> normals_;
  std::vector<Vector3D> positions_;
  std::vector<Vector2D> texCoords_;

  std::map<std::string, Material> materials_;
  std::vector<Object*> objects_;

  bool combined_;
};