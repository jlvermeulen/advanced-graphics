#include "renderViewer.h"
#include "ui_renderViewer.h"

#include <QFileDialog>
//--------------------------------------------------------------------------------
RenderViewer::RenderViewer(QWidget* parent)
  : QDialog(parent),
  ui(new Ui::RenderViewer)
{
  ui->setupUi(this);
}

//--------------------------------------------------------------------------------
RenderViewer::~RenderViewer()
{
  delete ui;
}

//--------------------------------------------------------------------------------
void RenderViewer::saveRender()
{
  QString filename = QFileDialog::getSaveFileName(parentWidget(), tr("Save rendered image"), "", tr("PNG Images (*.png)"));

  QFile file(filename);
  file.open(QIODevice::WriteOnly);

  ui->imageViewer->pixmap()->save(&file, "PNG");

  file.close();
}

//--------------------------------------------------------------------------------
void RenderViewer::setImage(QImage image)
{
  ui->imageViewer->setPixmap(QPixmap::fromImage(image));
}

//--------------------------------------------------------------------------------
void RenderViewer::setElaspedTime(int elapsedTime)
{
  QString time = QString("%1:%2.%3")
    .arg(elapsedTime / 60000, 2, 10, QChar('0'))
    .arg((elapsedTime % 60000) / 1000, 2, 10, QChar('0'))
    .arg(elapsedTime % 1000, 3, 10, QChar('0'));

  ui->elapsedTimeValue->setText(time);
}