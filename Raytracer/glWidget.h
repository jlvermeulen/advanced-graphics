#pragma once

#include "TriangleD.h"
#include "VertexD.h"

#include <QGLWidget>

typedef std::_Vector_iterator<std::_Vector_val<std::_Simple_types<std::string>>> IIterator;
typedef std::vector<Vector3D> Vector3DList;

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
  void parseLine(std::vector<std::string>& segments, Vector3DList& coords, Vector3DList& normals, Vector3DList& texCoords);
  void parseFace(IIterator& it, const IIterator& end, Vector3DList& coords, Vector3DList& normals, Vector3DList& texCoords);
  objType parseType(std::string& type);
  void parseLibrary(IIterator& it, const IIterator& end);
  void parseMaterial(IIterator& it, const IIterator& end);
  Vector3D parseNormal(IIterator& it) const;
  Vector3D parseSpace(IIterator& it) const;
  Vector3D parseTexCoords(IIterator& it) const;
  Vector3D parseVertex(IIterator& it) const;

  double parseDouble(const IIterator& iterator) const;
  int parseInteger(const IIterator& iterator) const;

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
  std::vector<TriangleD> triangles;
};
