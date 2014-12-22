#pragma once

#include <Triangle.h>
#include <Vertex.h>

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

  void initializeGL();
  void resizeGL(const int& w, const int& h);
  void paintGL();

  void mousePressEvent(QMouseEvent* event);
  void mouseMoveEvent(QMouseEvent* event);

  static void qNormalizeAngle(int& angle);

public slots:
  void setXRotation(int angle);
  void setYRotation(int angle);
  void setZRotation(int angle);

signals:
  void xRotationChanged(int angle);
  void yRotationChanged(int angle);
  void zRotationChanged(int angle);

private:
  int xRot;
  int yRot;
  int zRot;
  QPoint lastPos;
  std::vector<Triangle> triangles;
};
