#pragma once

#include <QMainWindow>

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

public:
  GLWidget* getGLWidget() const;

public slots:
  void openFileDialog();
  void openRenderDialog();

private:
  Ui::Window *ui;
};
