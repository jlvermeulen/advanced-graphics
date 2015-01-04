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

  bool getUseOctree() const;
  void setUseOctree(bool use);

  int getMinTriangles() const;
  void setMinTriangles(int minTriangles);

  int getMaxDepth() const;
  void setMaxDepth(int maxDepth);

private:
  Ui::RenderDialog* ui;
};