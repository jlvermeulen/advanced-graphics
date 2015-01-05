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
  QPoint getResolution() const { return QPoint(camera_.Width, camera_.Height);  }
  void setResolution(QPoint resolution){ camera_.Width = resolution.x(); camera_.Height = resolution.y(); }

  bool getUseOctree() const { return useOctree_;  }
  void setUseOctree(bool use) { useOctree_ = use; }

  int getMinTriangles() const { return minTriangles_; }
  void setMinTriangles(int minTriangles) { minTriangles_ = minTriangles; }

  int getMaxDepth() const { return maxDepth_; }
  void setMaxDepth(int maxDepth) { maxDepth_ = maxDepth; }

  void loadScene(QString& fileName);
  bool renderScene(uchar* imageData);

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
  struct Intersection
  {
    Intersection(const Vector3D& position_, const Vector3D& direction_, double time_, const Triangle& hit_) :
      hitPoint(position_ + time_ * direction_),
      hit(hit_)
    {
    }

    Vector3D hitPoint;
    Triangle hit;
  };

  struct Light
  {
    Light(const Vector3D& position_, const ColorD& color_) :
      position(position_),
      color(color_)
    {

    };

    Vector3D position;
    ColorD color;
  };
  void drawBoundingBoxes() const;
  void drawCameraRay() const;
  void drawLine(const Vector3D& v1, const Vector3D& v2) const;
  void drawModel();

  ColorD radiance(const Intersection& intersection, Ray ray, double refractiveIndex, int recursionDepth) const;
  ColorD traceRay(Ray ray, double refractiveIndex, int recursionDepth) const;

  ColorD calculateDiffuse(const Intersection& intersection) const;
  ColorD calculateReflection() const;
  ColorD calculateRefraction() const;

  // Rendering
  std::deque<Triangle> triangles;
  std::vector<Light> lights;
  Camera camera_;
  Ray debugRay_;
  bool boundingBoxVisible_;
  bool cameraRayVisible_;
  int recursionDepth;

  // Octree
  Octree octree;
  bool useOctree_;
  int minTriangles_;
  int maxDepth_;

  // Input
  float stepSize_;
  QPoint lastPos;
};
