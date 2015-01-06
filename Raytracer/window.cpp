#include "window.h"
#include "ui_window.h"

#include <QFileDialog>

#include "renderDialog.h"
#include "renderViewer.h"

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

GLWidget* Window::getGLWidget() const
{
  return ui->glWidget;
}

void Window::openFileDialog()
{
  QString fileName = QFileDialog::getOpenFileName(this, tr("Open Scene"), QString(), tr("Scene Files (*.obj)"));

  ui->glWidget->loadScene(fileName);
}

void Window::openRenderDialog()
{
  RenderDialog renderDialog(this);

  // Instantiate render dialog values
  renderDialog.setResolution(getGLWidget()->getResolution());
  renderDialog.setUseOctree(getGLWidget()->getUseOctree());
  renderDialog.setMinTriangles(getGLWidget()->getMinTriangles());
  renderDialog.setMaxDepth(getGLWidget()->getMaxDepth());
  renderDialog.setRayDistribution(getGLWidget()->getRayDistribution());
  renderDialog.setNumberOfRays(getGLWidget()->getNumberOfRays());
  renderDialog.setSigma(getGLWidget()->getSigma());

  if (renderDialog.exec())
  {
    QPoint resolution = renderDialog.getResolution();
    // Collect render dialog values
    getGLWidget()->setResolution(resolution);
    getGLWidget()->setUseOctree(renderDialog.getUseOctree());
    getGLWidget()->setMinTriangles(renderDialog.getMinTriangles());
    getGLWidget()->setMaxDepth(renderDialog.getMaxDepth());
    getGLWidget()->setRayDistribution(renderDialog.getRayDistribution());
    getGLWidget()->setNumberOfRays(renderDialog.getNumberOfRays());
    getGLWidget()->setSigma(renderDialog.getSigma());

    uchar* imageData = new uchar[resolution.x() * resolution.y() * 4];    // Width * Height * Color Channels
    int elapsedTime = getGLWidget()->renderScene(imageData);

    QImage image(imageData, resolution.x(), resolution.y(), QImage::Format_RGB32);

    RenderViewer renderViewer(this);

    // Instantiate render viewer image
    renderViewer.setImage(image);
    renderViewer.setElaspedTime(elapsedTime);

    renderViewer.exec();

    // Clean up
    delete imageData;
  }
}