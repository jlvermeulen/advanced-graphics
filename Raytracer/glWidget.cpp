#include "glWidget.h"

#include <QMouseEvent>

GLWidget::GLWidget(QWidget* parent)
  : QGLWidget(parent),
    xRot(0),
    yRot(0),
    zRot(0)
{

}

GLWidget::~GLWidget()
{

}

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

void GLWidget::resizeGL(const int& w, const int& h)
{
  glViewport(0, 0, (GLint) w, (GLint) h);

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glOrtho(-0.5, +0.5, -0.5, +0.5, 4.0, 15.0);
  glMatrixMode(GL_MODELVIEW);
}

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

void GLWidget::mousePressEvent(QMouseEvent* event)
{
  lastPos = event->pos();
}

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

static void qNormalizeAngle(int &angle)
{
    while (angle < 0)
        angle += 360 * 16;
    while (angle > 360 * 16)
        angle -= 360 * 16;
}

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
