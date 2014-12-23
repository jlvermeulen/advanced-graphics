#pragma once

#include <deque>
#include <Camera.h>
#include <Triangle.h>
#include <Vertex.h>
#include <Octree.h>

#include <QGLWidget>

class GLWidget : public QGLWidget
{
  Q_OBJECT

public:
  GLWidget(QWidget* parent);
  ~GLWidget();

public:
  void loadScene(QString& fileName);

protected:
  void glPerspective(double fovY, double aspect, double zNear, double zFar);

  void initializeGL();
  void resizeGL(int width, int height);
  void paintGL();

  void keyPressEvent(QKeyEvent* event);

  void mousePressEvent(QMouseEvent* event);
  void mouseMoveEvent(QMouseEvent* event);

private:
  Camera camera_;
  QPoint lastPos;
  std::deque<Triangle> triangles;
  Octree octree;
};
