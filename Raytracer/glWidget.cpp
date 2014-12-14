#include "glWidget.h"

#include "ObjReader.h"
#include "Vector3D.h"

#include <fstream>

#include <QMouseEvent>

//--------------------------------------------------------------------------------
GLWidget::GLWidget(QWidget* parent)
  : QGLWidget(parent),
    xRot(0),
    yRot(0),
    zRot(0)
{

}

//--------------------------------------------------------------------------------
GLWidget::~GLWidget()
{

}

//--------------------------------------------------------------------------------
void GLWidget::loadScene(QString& fileName)
{
  // Open file stream
  std::ifstream fin;
  fin.open(fileName.toUtf8().data());

  std::string line;

  std::vector<Vector3D> coords;
  std::vector<Vector3D> normals;
  std::vector<Vector3D> texCoords;

  while (std::getline(fin, line))
  {
    std::vector<std::string> segments = ObjReader::splitLine(line, ' ');

    parseLine(segments, coords, normals, texCoords);
  }
}

//--------------------------------------------------------------------------------
void GLWidget::parseLine(std::vector<std::string>& segments, Vector3DList& coords, Vector3DList& normals, Vector3DList& texCoords)
{
  IIterator it = segments.begin();
  IIterator end = segments.end();

  // No segments at all
  if (it == segments.end())
    return;

  switch (parseType(*it))
  {
    case objType::face:
      parseFace(++it, end, coords, normals, texCoords);
      break;

    case objType::library:
      parseLibrary(++it, end);
      break;

    case objType::material:
      parseMaterial(++it, end);
      break;

    case objType::normal:
      normals.push_back(parseNormal(++it));
      break;

    case objType::space:
      parseSpace(++it);
      break;

    case objType::texCoords:
      texCoords.push_back(parseTexCoords(++it));
      break;

    case objType::vertex:
      coords.push_back(parseVertex(++it));
      break;

    case objType::nil:
    default:
      // Don't do anything
      break;
  }
}

//--------------------------------------------------------------------------------
objType GLWidget::parseType(std::string& type)
{
  if (type == "v")
    return objType::vertex;
  else if (type == "vt")
    return objType::texCoords;
  else if (type == "vn")
    return objType::normal;
  else if (type == "vp")
    return objType::space;
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
void GLWidget::parseFace(IIterator& it, const IIterator& end, Vector3DList& coords, Vector3DList& normals, Vector3DList& texCoords)
{
  VertexD vertices[3];

  // Only parse triangles
  for (int i = 0; i < 3; ++i)
  {
    std::vector<std::string> indices = ObjReader::splitLine(*it, '/');

    IIterator vIt = indices.begin();

    VertexD vertex;
    vertex.Position = coords.at(parseInteger(vIt) - 1);
    vertex.UV = coords.at(parseInteger(++vIt) - 1);
    vertex.Normal = coords.at(parseInteger(++vIt) - 1);

    vertices[i] = vertex;
    ++it;
  }

  triangles.push_back(TriangleD(vertices));
}

//--------------------------------------------------------------------------------
void GLWidget::parseLibrary(IIterator& it, const IIterator& end)
{
  // TODO: mtllib import
}

//--------------------------------------------------------------------------------
void GLWidget::parseMaterial(IIterator& it, const IIterator& end)
{
  // TODO: usemtl
}

//--------------------------------------------------------------------------------
Vector3D GLWidget::parseNormal(IIterator& it) const
{
  Vector3D result;

  result.X = parseDouble(it);
  result.Y = parseDouble(++it);
  result.Z = parseDouble(++it);

  return result;
}

//--------------------------------------------------------------------------------
Vector3D GLWidget::parseSpace(IIterator& it) const
{ // Do we want to support this?
  Vector3D result;

  result.X = parseDouble(it);
  result.Y = parseDouble(++it);
  result.Z = parseDouble(++it);

  return result;
}

//--------------------------------------------------------------------------------
Vector3D GLWidget::parseTexCoords(IIterator& it) const
{
  Vector3D result;

  result.X = parseDouble(it);
  result.Y = parseDouble(++it);
  result.Z = parseDouble(++it);

  return result;
}

//--------------------------------------------------------------------------------
Vector3D GLWidget::parseVertex(IIterator& it) const
{
  Vector3D result;

  result.X = parseDouble(it);
  result.Y = parseDouble(++it);
  result.Z = parseDouble(++it);

  return result;
}

//--------------------------------------------------------------------------------
double GLWidget::parseDouble(const IIterator& it) const
{
  return atof((*it).c_str());
}

//--------------------------------------------------------------------------------
int GLWidget::parseInteger(const IIterator& it) const
{
  return atoi((*it).c_str());
}

//--------------------------------------------------------------------------------
void GLWidget::setXRotation(int angle)
{
  qNormalizeAngle(angle);

  if (angle != xRot)
  {
    xRot = angle;
    emit xRotationChanged(angle);
    updateGL();
  }
}

//--------------------------------------------------------------------------------
void GLWidget::setYRotation(int angle)
{
  qNormalizeAngle(angle);

  if (angle != yRot)
  {
    yRot = angle;
    emit xRotationChanged(angle);
    updateGL();
  }
}

//--------------------------------------------------------------------------------
void GLWidget::setZRotation(int angle)
{
  qNormalizeAngle(angle);

  if (angle != zRot)
  {
    zRot = angle;
    emit xRotationChanged(angle);
    updateGL();
  }
}

//--------------------------------------------------------------------------------
void GLWidget::initializeGL()
{
  qglClearColor(QColor::fromRgb(0, 0, 0, 255));

  //glClearColor(0.0, 0.0, 0.0, 0.0);

  glEnable(GL_DEPTH_TEST);
  glEnable(GL_CULL_FACE);
  glShadeModel(GL_SMOOTH);
  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  glEnable(GL_MULTISAMPLE);

  const GLfloat lightPosition[4] = { 0.5, 5.0, 7.0, 1.0 };
  glLightfv(GL_LIGHT0, GL_POSITION, lightPosition);
}

//--------------------------------------------------------------------------------
void GLWidget::resizeGL(const int& w, const int& h)
{
  glViewport(0, 0, (GLint) w, (GLint) h);

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glOrtho(-0.5, +0.5, -0.5, +0.5, 4.0, 15.0);
  glMatrixMode(GL_MODELVIEW);
}

//--------------------------------------------------------------------------------
void GLWidget::paintGL()
{
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  // Some random gl code for testing
  glLoadIdentity();
  glRotatef(xRot / 16.0, 1.0, 0.0, 0.0);
  glRotatef(yRot / 16.0, 0.0, 1.0, 0.0);
  glRotatef(zRot / 16.0, 0.0, 0.0, 1.0);

  // Draw triangle
  glBegin(GL_TRIANGLES);
  glColor3f(0.1, 0.2, 0.3);
  glVertex3f(-0.5, 0, 0);
  glVertex3f(+0.5, 0, 0);
  glVertex3f(0, +0.5, 0);
  glEnd();
}

//--------------------------------------------------------------------------------
void GLWidget::mousePressEvent(QMouseEvent* event)
{
  lastPos = event->pos();
}

//--------------------------------------------------------------------------------
void GLWidget::mouseMoveEvent(QMouseEvent* event)
{
  int dx = event->x() - lastPos.x();
  int dy = event->y() - lastPos.y();

  if (event->buttons() & Qt::LeftButton)
  {
    setXRotation(xRot + 8 * dy);
    setYRotation(yRot + 8 * dx);
  }
  else if (event->buttons() & Qt::RightButton)
  {
    setXRotation(xRot + 8 * dy);
    setZRotation(zRot + 8 * dx);
  }
}

//--------------------------------------------------------------------------------
void GLWidget::qNormalizeAngle(int& angle)
{
    while (angle < 0)
        angle += 360 * 16;
    while (angle > 360 * 16)
        angle -= 360 * 16;
}