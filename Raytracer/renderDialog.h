#pragma once

#include <QDialog>

namespace Ui {
  class RenderDialog;
}

class RenderDialog : public QDialog
{
  Q_OBJECT

public:
  explicit RenderDialog(QWidget* parent = 0);
  ~RenderDialog();

public:
  QPoint getResolution() const;
  void setResolution(QPoint resolution);

  int getMinTriangles() const;
  void setMinTriangles(int minTriangles);

  int getMaxDepth() const;
  void setMaxDepth(int maxDepth);

  int getNumberOfRays() const;
  void setNumberOfRays(int rays);

  double getSigma() const;
  void setSigma(double sigma);

private:
  Ui::RenderDialog* ui;
};