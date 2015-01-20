#include "renderDialog.h"
#include "ui_renderDialog.h"

//--------------------------------------------------------------------------------
RenderDialog::RenderDialog(QWidget* parent)
  : QDialog(parent),
    ui(new Ui::RenderDialog)
{
  ui->setupUi(this);

  // Fill resolution box
  ui->resolutionBox->addItem(QString("1980x1080"), QVariant(QPoint(1980, 1080)));
  ui->resolutionBox->addItem(QString("1280x720"), QVariant(QPoint(1280, 720)));
}

//--------------------------------------------------------------------------------
RenderDialog::~RenderDialog()
{
  delete ui;
}

//--------------------------------------------------------------------------------
QPoint RenderDialog::getResolution() const
{
  return ui->resolutionBox->currentData().toPoint();
}

//--------------------------------------------------------------------------------
void RenderDialog::setResolution(QPoint resolution)
{
  int index = ui->resolutionBox->findData(QVariant(resolution));
  ui->resolutionBox->setCurrentIndex(index);
}

//--------------------------------------------------------------------------------
int RenderDialog::getMinTriangles() const
{
  return ui->octreeMinTrianglesBox->value();
}

//--------------------------------------------------------------------------------
void RenderDialog::setMinTriangles(int minTriangles)
{
  ui->octreeMinTrianglesBox->setValue(minTriangles);
}

//--------------------------------------------------------------------------------
int RenderDialog::getMaxDepth() const
{
  return ui->octreeMaxDepthBox->value();
}

//--------------------------------------------------------------------------------
void RenderDialog::setMaxDepth(int maxDepth)
{
  ui->octreeMaxDepthBox->setValue(maxDepth);
}

//--------------------------------------------------------------------------------
int RenderDialog::getNumberOfRays() const
{
  return ui->numberOfRaysBox->value();
}

//--------------------------------------------------------------------------------
void RenderDialog::setNumberOfRays(int rays)
{
  ui->numberOfRaysBox->setValue(rays);
}

//--------------------------------------------------------------------------------
double RenderDialog::getSigma() const
{
  return ui->sigmaBox->value();
}

//--------------------------------------------------------------------------------
void RenderDialog::setSigma(double sigma)
{
  ui->sigmaBox->setValue(sigma);
}

//--------------------------------------------------------------------------------
bool RenderDialog::getUseDoF() const
{
  return ui->dofOnButton->isChecked();
}

//--------------------------------------------------------------------------------
void RenderDialog::setUseDoF(bool use)
{
  ui->dofOnButton->setChecked(use);
}

//--------------------------------------------------------------------------------
double RenderDialog::getFocalDepth() const
{
  return ui->focalDepthBox->value();
}

//--------------------------------------------------------------------------------
void RenderDialog::setFocalDepth(double distance)
{
  ui->focalDepthBox->setValue(distance);
}

//--------------------------------------------------------------------------------
double RenderDialog::getAperture() const
{
  return ui->apertureBox->value();
}

//--------------------------------------------------------------------------------
void RenderDialog::setAperture(double aperture)
{
  ui->apertureBox->setValue(aperture);
}