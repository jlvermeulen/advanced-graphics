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
  cameraControls = new CameraControls(this);
}

Window::~Window()
{
  delete cameraControls;
  delete ui;
}

GLWidget* Window::getGLWidget() const
{
  return ui->glWidget;
}

void Window::openCameraControls()
{
  cameraControls->setModal(false);

  cameraControls->setX(getGLWidget()->getX());
  cameraControls->setY(getGLWidget()->getY());
  cameraControls->setZ(getGLWidget()->getZ());

  connect(cameraControls->getXBox(), SIGNAL(valueChanged(double)),
          this, SLOT(setX(double)));

  connect(cameraControls->getYBox(), SIGNAL(valueChanged(double)),
          this, SLOT(setY(double)));
  
  connect(cameraControls->getZBox(), SIGNAL(valueChanged(double)),
          this, SLOT(setZ(double)));

  cameraControls->show();
  cameraControls->raise();
  cameraControls->activateWindow();
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
  renderDialog.setMinTriangles(getGLWidget()->getMinTriangles());
  renderDialog.setMaxDepth(getGLWidget()->getMaxDepth());
  renderDialog.setNumberOfRays(getGLWidget()->getNumberOfRays());
  renderDialog.setSigma(getGLWidget()->getSigma());
  renderDialog.setUseDoF(getGLWidget()->getUseDoF());
  renderDialog.setFocalDepth(getGLWidget()->getFocalDepth());
  renderDialog.setAperture(getGLWidget()->getAperture());

  if (renderDialog.exec())
  {
    QPoint resolution = renderDialog.getResolution();
    // Collect render dialog values
    getGLWidget()->setResolution(resolution);
    getGLWidget()->setMinTriangles(renderDialog.getMinTriangles());
    getGLWidget()->setMaxDepth(renderDialog.getMaxDepth());
    getGLWidget()->setNumberOfRays(renderDialog.getNumberOfRays());
    getGLWidget()->setSigma(renderDialog.getSigma());
    getGLWidget()->setUseDoF(renderDialog.getUseDoF());
    getGLWidget()->setFocalDistance(renderDialog.getFocalDepth());
    getGLWidget()->setAperture(renderDialog.getAperture());

    uchar* imageData = new uchar[resolution.x() * resolution.y() * 4];    // Width * Height * Color Channels
    int elapsedTime = getGLWidget()->renderScene(imageData);

    QImage image(imageData, resolution.x(), resolution.y(), QImage::Format_RGB32);

    RenderViewer renderViewer(this);

    // Instantiate render viewer image
    renderViewer.setImage(image);
    renderViewer.setElaspedTime(elapsedTime);

    renderViewer.exec();

    // Clean up
    delete [] imageData;
  }
}

void Window::setX(double x)
{
  getGLWidget()->setX(x);

  applyChanges();
}

void Window::setY(double y)
{
  getGLWidget()->setY(y);

  applyChanges();
}

void Window::setZ(double z)
{
  getGLWidget()->setZ(z);

  applyChanges();
}

void Window::applyChanges()
{
  activateWindow();
  cameraControls->activateWindow();
}