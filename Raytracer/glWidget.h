#pragma once

#include <QGLWidget>

class GLWidget : public QGLWidget
{
  Q_OBJECT

public:
  GLWidget(QWidget* parent);
  ~GLWidget();

public slots:
  void loadScene(QString& fileName);
  void setXRotation(int angle);
  void setYRotation(int angle);
  void setZRotation(int angle);

signals:
  void xRotationChanged(int angle);
  void yRotationChanged(int angle);
  void zRotationChanged(int angle);

protected:
  void initializeGL();
  void resizeGL(const int& w, const int& h);
  void paintGL();

  void mousePressEvent(QMouseEvent* event);
  void mouseMoveEvent(QMouseEvent* event);

private:
  int xRot;
  int yRot;
  int zRot;
  QPoint lastPos;
};
