#pragma once

#include <QDialog>
#include <QTime>


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
  void setElaspedTime(QTime time);

public slots:
  void saveRender();

private:
  Ui::RenderViewer* ui;
};