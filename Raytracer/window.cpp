#include "window.h"
#include "ui_window.h"

#include <QFileDialog>

Window::Window(QWidget *parent)
  : QMainWindow(parent),
    ui(new Ui::Window)
{
  ui->setupUi(this);
}

Window::~Window()
{
  delete ui;
}

void Window::openFileDialog()
{
  QString fileName = QFileDialog::getOpenFileName(this, tr("Open Scene"), QString(), tr("Scene Files (*.obj)"));

  ui->glWidget->loadScene(fileName);
}
