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
  ui->resolutionBox->addItem(QString("1280x720"), QVariant(QPoint(1280, 780)));

  // Fill data structure box
  ui->dataStructureBox->addItem(QString("List"), QVariant(false));
  ui->dataStructureBox->addItem(QString("Octree"), QVariant(true));
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
bool RenderDialog::getUseOctree() const
{
  return ui->dataStructureBox->currentData().toBool();
}

//--------------------------------------------------------------------------------
void RenderDialog::setUseOctree(bool use)
{
  int index = ui->dataStructureBox->findData(QVariant(use));
  ui->dataStructureBox->setCurrentIndex(index);
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