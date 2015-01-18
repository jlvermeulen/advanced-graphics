#include "cameraControls.h"
#include "ui_cameraControls.h"

//--------------------------------------------------------------------------------
CameraControls::CameraControls(QWidget* parent):
  QDialog(parent),
  ui(new Ui::CameraControls)
{
  ui->setupUi(this);
}

//--------------------------------------------------------------------------------
CameraControls::~CameraControls()
{
  delete ui;
}

//--------------------------------------------------------------------------------
QDoubleSpinBox* CameraControls::getXBox() const
{
  return ui->xBox;
}

//--------------------------------------------------------------------------------
QDoubleSpinBox* CameraControls::getYBox() const
{
  return ui->yBox;
}

//--------------------------------------------------------------------------------
QDoubleSpinBox* CameraControls::getZBox() const
{
  return ui->zBox;
}

//--------------------------------------------------------------------------------
double CameraControls::getX() const
{
  return ui->xBox->value();
}

//--------------------------------------------------------------------------------
void CameraControls::setX(double x)
{
  ui->xBox->setValue(x);
}

//--------------------------------------------------------------------------------
double CameraControls::getY() const
{
  return ui->yBox->value();
}

//--------------------------------------------------------------------------------
void CameraControls::setY(double y)
{
  ui->yBox->setValue(y);
}

//--------------------------------------------------------------------------------
double CameraControls::getZ() const
{
  return ui->zBox->value();
}

//--------------------------------------------------------------------------------
void CameraControls::setZ(double z)
{
  ui->zBox->setValue(z);
}

