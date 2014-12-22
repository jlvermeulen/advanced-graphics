#pragma once

#include <vector>
#include <Triangle.h>
#include <Vector3D.h>

typedef std::_String_const_iterator<std::_String_val<std::_Simple_types<char>>> CIIterator;
typedef std::_Vector_iterator<std::_Vector_val<std::_Simple_types<std::string>>> IIterator;
typedef std::vector<Vector3D> Vector3List;

enum objType
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

class ObjReader
{
public:
  ObjReader();
  ~ObjReader();

public:
  std::vector<Triangle> parseFile(const char* fileName);
  void parseLine(std::vector<std::string>& segments);

  void parseLibrary(IIterator& it, const IIterator& end);
  void parseMaterial(IIterator& it, const IIterator& end);
  void parseFace(IIterator& it);
  void parseNormal(IIterator& it);
  void parseTexCoords(IIterator& it);
  void parseVertex(IIterator& it);

  double parseDouble(const IIterator& iterator) const;
  int parseInteger(const IIterator& iterator) const;
  objType parseType(std::string& type);

  static std::vector<std::string> split(const std::string& text, char sep, bool multiple = true);

private:
  void reset();

private:
  std::vector<Vector3D> normals;
  std::vector<Vector3D> positions;
  std::vector<Vector3D> texCoords;
  std::vector<Triangle> triangles;

  bool combined_;
};