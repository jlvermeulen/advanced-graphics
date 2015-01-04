#define _USE_MATH_DEFINES

#include <glWidget.h>

#include <fstream>
#include <gl/GLU.h>
#include <math.h>
#include <ObjReader.h>
#include <QKeyEvent>
#include <QMouseEvent>
#include <QProgressDialog>
#include <Vector3D.h>
#include <queue>
#include <Intersections.h>

//--------------------------------------------------------------------------------
GLWidget::GLWidget(QWidget* parent)
  : QGLWidget(parent),
    boundingBoxVisible_(false),
    cameraRayVisible_(false),
    resolution_(QPoint(1980, 1080)),
    useOctree_(false),
    minTriangles_(0),
    maxDepth_(0),
    stepSize_(0.1f),
    lastPos()
{
}

//--------------------------------------------------------------------------------
GLWidget::~GLWidget()
{

}

//--------------------------------------------------------------------------------
void GLWidget::loadScene(QString& fileName)
{
  // Remove possible old data
  triangles.clear();

  ObjReader reader;
  triangles = reader.parseFile(fileName.toUtf8().data());

  if (useOctree_)
    octree = Octree(triangles, 100, 5);
}

//--------------------------------------------------------------------------------
bool GLWidget::renderScene(uchar* imageData) const
{
  // Show rendering progress
  //int numIterations = resolution_.x() * resolution_.y() * 3;  // Width * Height * Color Channels

  //QProgressDialog progressDialog("Rendering in progress...", "Cancel", 0, numIterations, parentWidget());
  //progressDialog.setWindowModality(Qt::WindowModal);

  QColor temp;
  Vector3D color = Vector3D();

  //int i = 0;

  for (int y = 0; y < resolution_.y(); ++y)
  {

    for (int x = 0; x < resolution_.x(); ++x)
    {
      int offset = (y * resolution_.x() + x) * 4;

      // For each color channel in reversed order (i.e. blue-green-red)
      for (int c = 0; c < 3; ++c)
      {
        //if (progressDialog.wasCanceled())
        //  return false;
        imageData[offset + c] = raytrace(x, y, c);

        //progressDialog.setValue(++i);
      }

      // Alpha channel
      imageData[offset + 3] = 255;
    }
  }

  //progressDialog.close();

  return true;
}

//--------------------------------------------------------------------------------
void GLWidget::setBoundingBoxVisible(bool visible)
{
  boundingBoxVisible_ = visible;
  updateGL();
}

//--------------------------------------------------------------------------------
void GLWidget::setCameraRayVisible(bool visible)
{
  cameraRayVisible_ = visible;
  debugRay_ = Ray(camera_.getEye(), camera_.getFocus(), ColorD(1.0, 1.0, 1.0));
  updateGL();
}

//--------------------------------------------------------------------------------
void GLWidget::initializeGL()
{
  //qglClearColor(QColor::fromRgb(0, 0, 0, 255));
  
  glClearColor(0.0, 0.0, 0.0, 0.0);

  glEnable(GL_DEPTH_TEST);
  glEnable(GL_CULL_FACE);
  glShadeModel(GL_SMOOTH);
  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  glEnable(GL_MULTISAMPLE);
  
  const GLfloat lightPosition[4] = { 0.5, 5.0, 7.0, 1.0 };
  glLightfv(GL_LIGHT0, GL_POSITION, lightPosition);
}

//--------------------------------------------------------------------------------
void GLWidget::resizeGL(int width, int height)
{
  glViewport(0, 0, width, height);

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glPerspective(40.0, (double) width / (double) height, 0.1, 20.0);

  glMatrixMode(GL_MODELVIEW);
}

//--------------------------------------------------------------------------------
void GLWidget::paintGL()
{
	if (triangles.empty())
		return;

  glEnable(GL_LIGHTING);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  Vector3D eye = camera_.getEye();
  Vector3D focus = camera_.getFocus();
  Vector3D viewPoint = eye + focus;
  Vector3D up = camera_.getUp();

  // Camera position
  glLoadIdentity();
  gluLookAt(eye.X, eye.Y, eye.Z,
            viewPoint.X, viewPoint.Y, viewPoint.Z,
            up.X, up.Y, up.Z);

  drawModel();

  // Turn off lights
  glDisable(GL_LIGHTING);

  if (useOctree_ && boundingBoxVisible_)
    drawBoundingBoxes();

  if (cameraRayVisible_)
    drawCameraRay();
}

//--------------------------------------------------------------------------------
void GLWidget::keyPressEvent(QKeyEvent* event)
{
  bool changed = false;

  if (event->key() == Qt::Key_W)        // Forward
  {
    changed = true;
    camera_.MoveForward(stepSize_);
  }
  else if (event->key() == Qt::Key_A)   // Left
  {
    changed = true;
    camera_.MoveRight(-stepSize_);
  }
  else if (event->key() == Qt::Key_S)   // Backward
  {
    changed = true;
    camera_.MoveForward(-stepSize_);
  }
  else if (event->key() == Qt::Key_D)   // Right
  {
    changed = true;
    camera_.MoveRight(stepSize_);
  }
  else if (event->key() == Qt::Key_X)   // Up
  {
    changed = true;
    camera_.MoveUpward(stepSize_);
  }
  else if (event->key() == Qt::Key_Z)   // Down
  {
    changed = true;
    camera_.MoveUpward(-stepSize_);
  }
  else
  {
    QGLWidget::keyPressEvent(event);
  }

  if (changed)
    updateGL();
}

//--------------------------------------------------------------------------------
void GLWidget::mousePressEvent(QMouseEvent* event)
{
  lastPos = event->pos();
}

//--------------------------------------------------------------------------------
void GLWidget::mouseMoveEvent(QMouseEvent* event)
{
  if (event->buttons() & Qt::LeftButton)
  {
    QPoint pos = event->pos();

    int dx = pos.x() - lastPos.x();
    int dy = pos.y() - lastPos.y();

    camera_.RotateY(0.1 * dx);         // Horizontal
    camera_.RotateX(0.1 * dy);         // Vertical

    updateGL();

    lastPos = pos;
  }
}

//--------------------------------------------------------------------------------
void GLWidget::glPerspective(double fovY, double aspect, double zNear, double zFar)
{
  double xMin, xMax, yMin, yMax;

  yMax = zNear * tan(fovY * M_PI / 360.0);
  yMin = -yMax;
  xMin = yMin * aspect;
  xMax = yMax * aspect;

  glFrustum(xMin, xMax, yMin, yMax, zNear, zFar);
}

//--------------------------------------------------------------------------------
void GLWidget::drawBoundingBoxes() const
{
  glEnableClientState(GL_VERTEX_ARRAY);
  //glBegin(GL_LINES);

  // Initialize queue
  std::queue<OctreeNode> q;
  q.push(octree.root);

  glColor3f(0, 1, 1);

  GLubyte indices[] = { // 24 indices
    0, 1, 0, 2,
    0, 4, 1, 3,
    1, 5, 2, 3,
    2, 6, 3, 7,
    4, 5, 4, 6,
    5, 7, 6, 7
  };

  while (!q.empty())
  {
    OctreeNode node = q.front();
    q.pop();

    BoundingBox bb = node.bb;
    //Vector3D verts[8];

    //for (int i = 0; i < 8; i++)
    //{
    //  double x = i & 4 ? bb.Halfsize.X : -bb.Halfsize.X;
    //  double y = i & 2 ? bb.Halfsize.Y : -bb.Halfsize.Y;
    //  double z = i & 1 ? bb.Halfsize.Z : -bb.Halfsize.Z;

    //  verts[i] = bb.Center + Vector3D(x, y, z);
    //}

    //drawLine(verts[0], verts[1]);
    //drawLine(verts[0], verts[2]);
    //drawLine(verts[0], verts[4]);
    //drawLine(verts[1], verts[3]);
    //drawLine(verts[1], verts[5]);
    //drawLine(verts[2], verts[3]);
    //drawLine(verts[2], verts[6]);
    //drawLine(verts[3], verts[7]);
    //drawLine(verts[4], verts[5]);
    //drawLine(verts[4], verts[6]);
    //drawLine(verts[5], verts[7]);
    //drawLine(verts[6], verts[7]);

    GLfloat vertices[24]; // 8 times 3 coordinates

    for (int i = 0; i < 8; i++)
    {
      float x = i & 4 ? bb.Halfsize.X : -bb.Halfsize.X;
      float y = i & 2 ? bb.Halfsize.Y : -bb.Halfsize.Y;
      float z = i & 1 ? bb.Halfsize.Z : -bb.Halfsize.Z;

      vertices[i * 3] = bb.Center.X + x;
      vertices[i * 3 + 1] = bb.Center.Y + y;
      vertices[i * 3 + 2] = bb.Center.Z + z;
    }

    glVertexPointer(3, GL_FLOAT, 0, vertices);

    glDrawElements(GL_LINES, 24, GL_UNSIGNED_BYTE, indices);

    for (OctreeNode& n : node.children)
      q.push(n);
  }

  glDisableClientState(GL_VERTEX_ARRAY);
  //glEnd();
}

//--------------------------------------------------------------------------------
void GLWidget::drawCameraRay() const
{
  glBegin(GL_LINES);

  Triangle tri;
  double time;
  bool hit = octree.Query(debugRay_, tri, time);

  if (hit)
  {
    Vector3D point = debugRay_.Origin + time * debugRay_.Direction;
    glColor3f(1, 0, 0);
    double eps = 0.0025;
    for (int i = 0; i < 8; i++)
    {
      double x = i & 1 ? eps : -eps;
      double y = i & 2 ? eps : -eps;
      double z = i & 4 ? eps : -eps;
      drawLine(point, Vector3D(point.X + x, point.Y + y, point.Z + z));
    }
    glColor3f(0, 1, 0);
  }
  else
  {
    glColor3f(1, 0, 0);
  }

  drawLine(debugRay_.Origin, debugRay_.Origin + 10 * debugRay_.Direction);

  glEnd();
}

//--------------------------------------------------------------------------------
void GLWidget::drawLine(const Vector3D& v1, const Vector3D& v2) const
{
  glVertex3f(v1.X, v1.Y, v1.Z);
  glVertex3f(v2.X, v2.Y, v2.Z);
}

//--------------------------------------------------------------------------------
void GLWidget::drawModel()
{
  glBegin(GL_TRIANGLES);

  for (Triangle& triangle : triangles)
  {
    for (Vertex& vertex : triangle.Vertices)
    {
      glColor3f(vertex.Color.R, vertex.Color.G, vertex.Color.B);
      glNormal3f(vertex.Normal.X, vertex.Normal.Y, vertex.Normal.Z);
      glVertex3f(vertex.Position.X, vertex.Position.Y, vertex.Position.Z);
    }
  }

  glEnd();
}

//--------------------------------------------------------------------------------
uchar GLWidget::raytrace(int x, int y, int channel) const
{


  return 0; // between 0 and 255
}