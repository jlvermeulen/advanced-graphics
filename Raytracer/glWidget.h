#pragma once

#include <deque>
#include <QGLWidget>

#include "Camera.h"
#include "Triangle.h"
#include "Vertex.h"
#include "Octree.h"


class GLWidget : public QGLWidget
{
  Q_OBJECT

public:
  GLWidget(QWidget* parent);
  ~GLWidget();

public:
  QPoint getResolution() const { return resolution_;  }
  void setResolution(QPoint resolution) { resolution_ = resolution; }

  bool getUseOctree() const { return useOctree_;  }
  void setUseOctree(bool use) { useOctree_ = use; }

  int getMinTriangles() const { return minTriangles_; }
  void setMinTriangles(int minTriangles) { minTriangles_ = minTriangles; }

  int getMaxDepth() const { return maxDepth_; }
  void setMaxDepth(int maxDepth) { maxDepth_ = maxDepth; }

  void loadScene(QString& fileName);
  bool renderScene(uchar* imageData) const;

public slots:
  void setBoundingBoxVisible(bool visible);
  void setCameraRayVisible(bool visible);

protected:
  void glPerspective(double fovY, double aspect, double zNear, double zFar);

  void initializeGL();
  void resizeGL(int width, int height);
  void paintGL();

  void keyPressEvent(QKeyEvent* event);

  void mousePressEvent(QMouseEvent* event);
  void mouseMoveEvent(QMouseEvent* event);

private:
  void drawBoundingBoxes() const;
  void drawCameraRay() const;
  void drawLine(const Vector3D& v1, const Vector3D& v2) const;
  void drawModel();

  uchar raytrace(int x, int y, int channel) const;

  // Rendering
  std::deque<Triangle> triangles;
  Camera camera_;
  Ray debugRay_;
  bool boundingBoxVisible_;
  bool cameraRayVisible_;
  QPoint resolution_;

  // Octree
  Octree octree;
  bool useOctree_;
  int minTriangles_;
  int maxDepth_;

  // Input
  float stepSize_;
  QPoint lastPos;
};
