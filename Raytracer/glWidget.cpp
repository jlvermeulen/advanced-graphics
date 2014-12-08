#include "glWidget.h"

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

  while (std::getline(fin, line))
    parseLine(line);
}

//--------------------------------------------------------------------------------
void GLWidget::parseLine(std::string& line)
{
  IIterator begin = line.begin();
  IIterator iterator = line.begin();
  IIterator end = line.end();

  // Read until first white-space or end-of-file
  while (iterator != line.end() && *iterator != ' ')
    ++iterator;

  // Corner-case where iterator == begin is not handled!
  std::string type(begin, iterator - 1);

  switch (parseType(type))
  {
    case objType::face:
      parseFace(++iterator, end);
      break;

    case objType::library:
      parseLibrary(++iterator, end);
      break;

    case objType::material:
      parseMaterial(++iterator, end);
      break;

    case objType::normal:
      parseNormal(++iterator, end);
      break;

    case objType::space:
      parseSpace(++iterator, end);
      break;

    case objType::texCoords:
      parseTexCoords(++iterator, end);
      break;

    case objType::vertex:
      parseVertex(++iterator, end);
      break;

    case objType::nil:
    default:
      // Don't do anything
      break;
  }

  //size_t length = line.length();

  //std::vector<std::string> args;
  //std::string arg;

  //size_t startPos = 0;
  //size_t endPos;

  //do 
  //{
  //  // Find end of argument
  //  endPos = line.find(' ', startPos);

  //  // Retrieve argument
  //  if (endPos == std::string::npos)
  //    arg = line.substr(startPos);
  //  else
  //    arg = line.substr(startPos, endPos - startPos);

  //  // Add argument
  //  args.push_back(arg);

  //  // Move to next character
  //  if (endPos < length)
  //    startPos = endPos + 1;
  //  else
  //    startPos = std::string::npos;
  //} while (startPos != std::string::npos);
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
void GLWidget::parseFace(IIterator& iterator, const IIterator& end)
{
  std::vector<int> vertices;
  std::vector<int> texCoords;
  std::vector<int> normals;

  while (false) // TODO: Read only until next letter
  {
    // TODO: make parts optional?
    vertices.push_back(parseInteger(iterator, end));
    texCoords.push_back(parseInteger(iterator, end));
    normals.push_back(parseInteger(iterator, end));
  }
}

//--------------------------------------------------------------------------------
void GLWidget::parseLibrary(IIterator& iterator, const IIterator& end)
{
  // TODO: mtllib import
}

//--------------------------------------------------------------------------------
void GLWidget::parseMaterial(IIterator& iterator, const IIterator& end)
{
  // TODO: usemtl
}

//--------------------------------------------------------------------------------
void GLWidget::parseNormal(IIterator& iterator, const IIterator& end)
{
  double x = parseDouble(iterator, end);
  double y = parseDouble(iterator, end);
  double z = parseDouble(iterator, end);
}

//--------------------------------------------------------------------------------
void GLWidget::parseSpace(IIterator& iterator, const IIterator& end)
{ // Do we want to support this?
  double u = parseDouble(iterator, end);
  double v = parseDouble(iterator, end);
  // double w = parseDouble(iterator, end);
  double w = 1.0;
}

//--------------------------------------------------------------------------------
void GLWidget::parseTexCoords(IIterator& iterator, const IIterator& end)
{
  double u = parseDouble(iterator, end);
  double v = parseDouble(iterator, end);
  // double w = parseDouble(iterator, end);
  double w = 1.0;
}

//--------------------------------------------------------------------------------
void GLWidget::parseVertex(IIterator& iterator, const IIterator& end)
{
  double x = parseDouble(iterator, end);
  double y = parseDouble(iterator, end);
  double z = parseDouble(iterator, end);
  // double w = parseDouble(iterator, end);
  double w = 1.0;
}

//--------------------------------------------------------------------------------
double GLWidget::parseDouble(IIterator& iterator, const IIterator& end)
{
  IIterator begin(iterator);

  // Read until first white-space or end-of-file
  while (iterator != end && *iterator != ' ')
    ++iterator;

  // Corner-case, where equal
  std::string value(begin, iterator - 1);

  return atof(value.c_str());
}

//--------------------------------------------------------------------------------
int GLWidget::parseInteger(IIterator& iterator, const IIterator& end)
{
  IIterator begin(iterator);

  // Read until first white-space or end-of-file
  while (iterator != end && *iterator != ' ' && *iterator != '/')
    ++iterator;

  // Corner-case, where equal
  std::string value(begin, iterator - 1);

  return atoi(value.c_str());
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