#pragma once

#include <QDialog>
#include <QDoubleSpinBox>

namespace Ui {
  class CameraControls;
}

class CameraControls : public QDialog
{
  Q_OBJECT

public:
  explicit CameraControls(QWidget* parent = 0);
  ~CameraControls();

public:
  QDoubleSpinBox* getXBox() const;
  QDoubleSpinBox* getYBox() const;
  QDoubleSpinBox* getZBox() const;

  double getX() const;
  void setX(double x);

  double getY() const;
  void setY(double y);

  double getZ() const;
  void setZ(double z);

private:
  Ui::CameraControls* ui;
};