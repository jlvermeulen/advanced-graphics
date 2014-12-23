#include <ObjReader.h>

#include <fstream>
#include <string>

//--------------------------------------------------------------------------------
ObjReader::ObjReader() :
  combined_(false)
{

}

//--------------------------------------------------------------------------------
ObjReader::~ObjReader()
{

}

//--------------------------------------------------------------------------------
std::deque<Triangle> ObjReader::parseFile(const char* fileName)
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

  return triangles;
}

//--------------------------------------------------------------------------------
void ObjReader::parseLine(std::vector<std::string>& segments)
{
  IIterator it = segments.begin();
  IIterator end = segments.end();

  // No segments at all
  if (it == segments.end())
    return;

  switch (parseType(*it))
  {
    case objType::face:
      // Last step of parsing an object
      parseFace(++it);
      combined_ = true;
      break;

    case objType::library:
      parseLibrary(++it, end);
      break;

    case objType::material:
      parseMaterial(++it, end);
      break;

    case objType::normal:
      parseNormal(++it);
      break;

    case objType::texCoords:
      parseTexCoords(++it);
      break;

    case objType::vertex:
      parseVertex(++it);
      break;

    case objType::nil:
    default:
    {
      // Don't do anything
      break;
    }
  }
}

//--------------------------------------------------------------------------------
void ObjReader::parseFace(IIterator& it)
{
  Vertex vertices[3];

  // Only parse triangles
  for (int i = 0; i < 3; ++i)
  {
    std::vector<std::string> indices = ObjReader::split(*it, '/');

    IIterator vIt = indices.begin();

    Vertex vertex;
    vertex.Position = positions.at(parseInteger(vIt) - 1);

    if (indices.size() == 3)
      vertex.UV = texCoords.at(parseInteger(++vIt) - 1);

    vertex.Normal = normals.at(parseInteger(++vIt) - 1);

    vertices[i] = vertex;
    ++it;
  }

  // Add triangle
  triangles.push_back(Triangle(vertices));
}

//--------------------------------------------------------------------------------
void ObjReader::parseLibrary(IIterator& it, const IIterator& end)
{
  // TODO: mtllib import
}

//--------------------------------------------------------------------------------
void ObjReader::parseMaterial(IIterator& it, const IIterator& end)
{
  // TODO: usemtl
}

//--------------------------------------------------------------------------------
void ObjReader::parseNormal(IIterator& it)
{
  Vector3D normal;

  normal.X = parseDouble(it);
  normal.Y = parseDouble(++it);
  normal.Z = parseDouble(++it);

  normals.push_back(normal);
}

//--------------------------------------------------------------------------------
void ObjReader::parseTexCoords(IIterator& it)
{
  Vector3D texCoord;

  texCoord.X = parseDouble(it);
  texCoord.Y = parseDouble(++it);
  texCoord.Z = parseDouble(++it);

  texCoords.push_back(texCoord);
}

//--------------------------------------------------------------------------------
void ObjReader::parseVertex(IIterator& it)
{
  Vector3D position;

  position.X = parseDouble(it);
  position.Y = parseDouble(++it);
  position.Z = parseDouble(++it);

  positions.push_back(position);
}

//--------------------------------------------------------------------------------
double ObjReader::parseDouble(const IIterator& it) const
{
  return atof((*it).c_str());
}

//--------------------------------------------------------------------------------
int ObjReader::parseInteger(const IIterator& it) const
{
  return atoi((*it).c_str());
}

//--------------------------------------------------------------------------------
objType ObjReader::parseType(std::string& type)
{
  if (type == "v")
    return objType::vertex;
  else if (type == "vt")
    return objType::texCoords;
  else if (type == "vn")
    return objType::normal;
  else if (type == "f")
    return objType::face;
  else if (type == "mtllib")
    return objType::library;
  else if (type == "usemtl")
    return objType::material;
  else
    return objType::nil;
}

//--------------------------------------------------------------------------------
void ObjReader::reset()
{
  normals.clear();
  positions.clear();
  texCoords.clear();

  combined_ = false;
}

//--------------------------------------------------------------------------------
std::vector<std::string> ObjReader::split(const std::string& text, char sep, bool multiple)
{
  std::vector<std::string> result;
  std::string part;
  bool set = false;
  bool prev = false;

  CIIterator it = text.cbegin();

  while (it != text.cend())
  {
    if (*it == sep)
    {
      if (set || (!multiple && prev))
      {
        result.push_back(part);
        part = "";

        set = false;
        prev = false;
      }
      else
      {
        prev = true;
      }
    }
    else
    {
      set = true;
      part += *it;
    }

    ++it;
  }

  // Add the remaining part
  if (part.size() > 0)
    result.push_back(part);

  return result;
}