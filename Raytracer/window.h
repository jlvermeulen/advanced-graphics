#pragma once

#include <QMainWindow>

namespace Ui {
  class Window;
}

class Window : public QMainWindow
{
  Q_OBJECT

public:
  explicit Window(QWidget *parent = 0);
  ~Window();

public slots:
  void openFileDialog();

private:
  Ui::Window *ui;
};
