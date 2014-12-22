#define _USE_MATH_DEFINES

#include <glWidget.h>

#include <fstream>
#include <math.h>
#include <ObjReader.h>
#include <QKeyEvent>
#include <QMouseEvent>
#include <Vector3D.h>

//--------------------------------------------------------------------------------
GLWidget::GLWidget(QWidget* parent)
  : QGLWidget(parent),
    lastPos()
{
}

//--------------------------------------------------------------------------------
GLWidget::~GLWidget()
{

}

//--------------------------------------------------------------------------------
void GLWidget::loadScene(QString& fileName)
{
  // Remove possible old data
  triangles.clear();
  triangles.reserve(0);

  ObjReader reader;
  triangles = reader.parseFile(fileName.toUtf8().data());
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
  glPerspective(40.0, (double) w / (double) h, 0.5, 20.0);
  glMatrixMode(GL_MODELVIEW);
}

//--------------------------------------------------------------------------------
void GLWidget::paintGL()
{
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  Vector3D position = camera_.getPosition();
  Vector3D rotation = camera_.getRotation();

  // Camera position
  glLoadIdentity();
  glRotatef(-rotation.X, 1.0, 0.0, 0.0);
  glRotatef(-rotation.Y, 0.0, 1.0, 0.0);
  glRotatef(-rotation.Z, 0.0, 0.0, 1.0);
  glTranslatef(-position.X, -position.Y, -position.Z);

  // Draw triangles
  glBegin(GL_TRIANGLES);

  for (Triangle& triangle : triangles)
  {
    for (Vertex& vertex : triangle.Vertices)
    {
      glColor3f(vertex.Color.R, vertex.Color.G, vertex.Color.B);
      glVertex3f(vertex.Position.X, vertex.Position.Y, vertex.Position.Z);
    }
  }

  glEnd();
}

//--------------------------------------------------------------------------------
void GLWidget::keyPressEvent(QKeyEvent* event)
{
  float step = 0.05;
  bool changed = false;

  if (event->key() == Qt::Key_W)        // Forward
  {
    changed = true;
    camera_.MoveForward(step);
  }
  else if (event->key() == Qt::Key_A)   // Left
  {
    changed = true;
    camera_.MoveSideways(-step);
  }
  else if (event->key() == Qt::Key_S)   // Backward
  {
    changed = true;
    camera_.MoveForward(-step);
  }
  else if (event->key() == Qt::Key_D)   // Right
  {
    changed = true;
    camera_.MoveSideways(step);
  }
  else if (event->key() == Qt::Key_Z)   // Up
  {
    changed = true;
    camera_.Move(Vector3D::Up);
  }
  else if (event->key() == Qt::Key_X)   // Down
  {
    changed = true;
    camera_.Move(-Vector3D::Up);
  }

  if (changed)
    updateGL();
}

//--------------------------------------------------------------------------------
void GLWidget::mouseMoveEvent(QMouseEvent* event)
{
  if (event->buttons() & Qt::LeftButton)
  {
    QPoint pos = event->pos();

    int dx = pos.x() - lastPos.x();
    int dy = pos.y() - lastPos.y();

    camera_.RotateX(0.1 * dx);         // Horizontal
    camera_.RotateY(0.1 * dy);         // Vertical

    updateGL();

    lastPos = pos;
  }
}



//--------------------------------------------------------------------------------
void GLWidget::glPerspective(double fovY, double aspect, double zNear, double zFar)
{
  double xMin, xMax, yMin, yMax;

  yMax = zNear * tan(fovY * M_PI / 360.0);
  yMin = -yMax;
  xMin = yMin * aspect;
  xMax = yMax * aspect;

  glFrustum(xMin, xMax, yMin, yMax, zNear, zFar);
}