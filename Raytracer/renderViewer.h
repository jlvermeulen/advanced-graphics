#pragma once

#include <QDialog>

namespace Ui {
  class RenderViewer;
}

class RenderViewer : public QDialog
{
  Q_OBJECT

public:
  explicit RenderViewer(QWidget* parent = 0);
  ~RenderViewer();

public:
  void setImage(QImage image);

public slots:
  void saveRender();

private:
  Ui::RenderViewer* ui;
};