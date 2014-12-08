#pragma once

#include <QGLWidget>

typedef std::_String_iterator<std::_String_val<std::_Simple_types<char>>> IIterator;

enum objType
{
  face,
  normal,
  space,
  texCoords,
  vertex,
  library,
  material,
  nil
};

class GLWidget : public QGLWidget
{
  Q_OBJECT

public:
  GLWidget(QWidget* parent);
  ~GLWidget();

public:
  void loadScene(QString& fileName);

protected:
  void parseLine(std::string& line);
  objType parseType(std::string& type);

  void parseFace(IIterator& iterator, const IIterator& end);
  void parseLibrary(IIterator& iterator, const IIterator& end);
  void parseMaterial(IIterator& iterator, const IIterator& end);
  void parseNormal(IIterator& iterator, const IIterator& end);
  void parseSpace(IIterator& iterator, const IIterator& end);
  void parseTexCoords(IIterator& iterator, const IIterator& end);
  void parseVertex(IIterator& iterator, const IIterator& end);

  double parseDouble(IIterator& iterator, const IIterator& end);
  int parseInteger(IIterator& iterator, const IIterator& end);

  void initializeGL();
  void resizeGL(const int& w, const int& h);
  void paintGL();

  void mousePressEvent(QMouseEvent* event);
  void mouseMoveEvent(QMouseEvent* event);

  static void qNormalizeAngle(int& angle);

public slots:
  void setXRotation(int angle);
  void setYRotation(int angle);
  void setZRotation(int angle);

signals:
  void xRotationChanged(int angle);
  void yRotationChanged(int angle);
  void zRotationChanged(int angle);

private:
  int xRot;
  int yRot;
  int zRot;
  QPoint lastPos;
};
