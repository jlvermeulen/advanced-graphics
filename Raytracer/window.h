#pragma once

#include <QMainWindow>

#include "cameraControls.h"
#include "glWidget.h"

namespace Ui {
  class Window;
}

class Window : public QMainWindow
{
  Q_OBJECT

public:
  explicit Window(QWidget *parent = 0);
  ~Window();

  GLWidget* getGLWidget() const;

public slots:
  void openCameraControls();
  void openFileDialog();
  void openRenderDialog();

  void setX(double x);
  void setY(double y);
  void setZ(double z);

private:
  void applyChanges();

  Ui::Window *ui;
  CameraControls* cameraControls;;
};
