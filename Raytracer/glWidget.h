#pragma once

#include <deque>
#include <QGLWidget>
#include <QTime>

#include "Camera.h"
#include "Triangle.h"
#include "Vertex.h"
#include "Octree.h"
#include "Scene.h"


class GLWidget : public QGLWidget
{
  Q_OBJECT

public:
  GLWidget(QWidget* parent);
  ~GLWidget();

public:
  QPoint getResolution() const { return QPoint(scene.camera.Width, scene.camera.Height);  }
  void setResolution(QPoint resolution){ scene.camera.Width = resolution.x(); scene.camera.Height = resolution.y(); }

  bool getUseOctree() const { return useOctree_;  }
  void setUseOctree(bool use) { useOctree_ = use; }

  int getMinTriangles() const { return minTriangles_; }
  void setMinTriangles(int minTriangles) { minTriangles_ = minTriangles; }

  int getMaxDepth() const { return maxDepth_; }
  void setMaxDepth(int maxDepth) { maxDepth_ = maxDepth; }

  int getRayDistribution() const { return rayDistribution_; }
  void setRayDistribution(int distribution) { rayDistribution_ = static_cast<RayDistributionType>(distribution); }

  int getNumberOfRays() const { return numberOfRays_; }
  void setNumberOfRays(int rays) { numberOfRays_ = rays; }

  double getSigma() const { return sigma_; }
  void setSigma(double sigma) { sigma_ = sigma; }

  void loadScene(QString& fileName);
  int renderScene(uchar* imageData);

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

  // Rendering
  Scene scene;
  Ray debugRay_;
  bool boundingBoxVisible_;
  bool cameraRayVisible_;
  int numberOfRays_;
  RayDistributionType rayDistribution_;
  double sigma_;

  // Octree
  bool useOctree_;
  int minTriangles_;
  int maxDepth_;

  // Input
  float stepSize_;
  QPoint lastPos;
};
