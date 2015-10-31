#include "ObjReader.h"

#include "MtlReader.h"

#include <fstream>
#include <string>

//--------------------------------------------------------------------------------
ObjReader::ObjReader() :
  Reader(),
  combined_(true)
{

}

//--------------------------------------------------------------------------------
ObjReader::~ObjReader()
{

}

//--------------------------------------------------------------------------------
std::vector<Object*> ObjReader::parseFile(const char* fileName)
{
  // Out with the old
  reset();
  
  std::string path = std::string(fileName);
  path = path.substr(0, path.find_last_of('\\/') + 1);

  // Open file stream
  std::ifstream fin;
  fin.open(fileName);

  std::string line;

  while (std::getline(fin, line))
  {
    std::vector<std::string> segments = split(line, ' ');

    parseLine(path, segments);
  }

  return objects_;
}

//--------------------------------------------------------------------------------
void ObjReader::parseLine(const std::string& path, const std::vector<std::string>& segments)
{
  CVSIterator it = segments.cbegin();

  // No segments at all
  if (it == segments.cend())
    return;

  switch (parseType(*it))
  {
    case ObjType::face:
    {
      // Last step of parsing an object
      parseFace(++it);
      combined_ = true;
      break;
    }

    case ObjType::library:
      parseLibrary(path, ++it);
      break;

    case ObjType::material:
      parseMaterial(++it);
      break;

    case ObjType::normal:
      parseNormal(++it);
      break;

    case ObjType::texCoords:
      parseTexCoords(++it);
      break;

    case ObjType::vertex:
      // Do we need a new object?
      if (combined_)
      {
        // Create empty object
        objects_.push_back(new Object());

        combined_ = false;
      }

      parseVertex(++it);
      break;

    case ObjType::nil:
    default:
      // Don't do anything
      break;
  }
}

//--------------------------------------------------------------------------------
void ObjReader::parseFace(CVSIterator& it)
{
  Vertex vertices[3];

  // Only parse triangles
  for (int i = 0; i < 3; ++i)
  {
    std::vector<std::string> indices = split(*it, '/');

    CVSIterator vIt = indices.cbegin();

    Vertex vertex;
    vertex.Position = positions_.at(parseInteger(vIt) - 1);

    if (indices.size() == 3)
      vertex.UV = texCoords_.at(parseInteger(++vIt) - 1);

    vertex.Normal = normals_.at(parseInteger(++vIt) - 1);

    vertices[i] = vertex;
    ++it;
  }

  // Add triangle
  objects_.back()->triangles.push_back(new Triangle(vertices));
}

//--------------------------------------------------------------------------------
void ObjReader::parseLibrary(const std::string& path, CVSIterator& it)
{
  std::string fileName = *it;

  MtlReader mtlReader;

  materials_ = mtlReader.parseFile((path + fileName).c_str());
}

//--------------------------------------------------------------------------------
void ObjReader::parseMaterial(CVSIterator& it)
{
  std::string name = *it;

  objects_.back()->material = materials_.at(name);
}

//--------------------------------------------------------------------------------
void ObjReader::parseNormal(CVSIterator& it)
{
  Vector3D normal;

  normal.X = parseDouble(it);
  normal.Y = parseDouble(++it);
  normal.Z = parseDouble(++it);

  normals_.push_back(normal);
}

//--------------------------------------------------------------------------------
void ObjReader::parseTexCoords(CVSIterator& it)
{
  Vector2D texCoord;

  texCoord.X = parseDouble(it);
  texCoord.Y = parseDouble(++it);

  texCoords_.push_back(texCoord);
}

//--------------------------------------------------------------------------------
void ObjReader::parseVertex(CVSIterator& it)
{
  Vector3D position;

  position.X = parseDouble(it);
  position.Y = parseDouble(++it);
  position.Z = parseDouble(++it);

  positions_.push_back(position);
}

//--------------------------------------------------------------------------------
ObjType::ObjType ObjReader::parseType(const std::string& type) const
{
  if (type == "v")
    return ObjType::vertex;
  else if (type == "vt")
    return ObjType::texCoords;
  else if (type == "vn")
    return ObjType::normal;
  else if (type == "f")
    return ObjType::face;
  else if (type == "mtllib")
    return ObjType::library;
  else if (type == "usemtl")
    return ObjType::material;
  else
    return ObjType::nil;
}

//--------------------------------------------------------------------------------
void ObjReader::reset()
{
  normals_.clear();
  positions_.clear();
  texCoords_.clear();
}