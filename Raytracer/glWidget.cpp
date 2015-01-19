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

//--------------------------------------------------------------------------------
GLWidget::GLWidget(QWidget* parent)
  : QGLWidget(parent),
    boundingBoxVisible_(false),
    cameraRayVisible_(false),
    numberOfRays_(1),
    sigma_(0.1),
    minTriangles_(10),
    maxDepth_(10),
    stepSize_(0.1f),
    lastPos(),
    scene()
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
  scene.objects.clear();

  ObjReader reader;

  scene.objects = reader.parseFile(fileName.toUtf8().data());

  // Add lights
  //scene.lights.push_back(Light(Vector3D(-3.0, -5.0, -4.0), ColorD(15.0, 15.0, 15.0)));
  //scene.lights.push_back(Light(Vector3D(3.0, 5.0, 4.0), ColorD(15.0, 15.0, 15.0)));
}

//--------------------------------------------------------------------------------
int GLWidget::renderScene(uchar* imageData)
{
  QTime timer;
  timer.start();

  scene.Render(imageData, minTriangles_, maxDepth_, numberOfRays_, sigma_);

  return timer.elapsed();

  // Progress dialog slows it down enormously, so left it out..
  // Show rendering progress
  //int numIterations = scene.camera.Width * scene.camera.Height;

  //QProgressDialog progressDialog("Rendering in progress...", "Cancel", 0, numIterations, parentWidget());
  //progressDialog.setWindowModality(Qt::WindowModal);
  //progressDialog.close();
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
  debugRay_ = Ray(scene.camera.Eye(), scene.camera.Focus());
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
  //const GLfloat lightDiffuse[4] = { 1.0, 1.0, 1.0, 1.0 };
  glLightfv(GL_LIGHT0, GL_POSITION, lightPosition);
  //glLightfv(GL_LIGHT0, GL_DIFFUSE, lightDiffuse);
}

//--------------------------------------------------------------------------------
void GLWidget::resizeGL(int width, int height)
{
  glViewport(0, 0, width, height);

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glPerspective(scene.camera.FovY(), (double) width / (double) height, scene.camera.ZNear(), scene.camera.ZFar());

  glMatrixMode(GL_MODELVIEW);
}

//--------------------------------------------------------------------------------
void GLWidget::paintGL()
{
  glEnable(GL_LIGHTING);
  glEnable(GL_COLOR_MATERIAL);
  glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  Vector3D viewPoint = scene.camera.Eye() + scene.camera.Focus();

  // Camera position
  glLoadIdentity();
  gluLookAt(scene.camera.Eye().X, scene.camera.Eye().Y, scene.camera.Eye().Z,
            viewPoint.X, viewPoint.Y, viewPoint.Z,
            scene.camera.Up().X, scene.camera.Up().Y, scene.camera.Up().Z);

  drawModel();

  // Turn off lights
  glDisable(GL_LIGHTING);

  if (boundingBoxVisible_)
    drawBoundingBoxes();

  if (cameraRayVisible_)
    drawCameraRay();
}

//--------------------------------------------------------------------------------
void GLWidget::keyPressEvent(QKeyEvent* event)
{
  bool changed = true;

  if (event->key() == Qt::Key_W)        // Forward
    scene.camera.MoveForward(stepSize_);
  else if (event->key() == Qt::Key_A)   // Left
    scene.camera.MoveRight(-stepSize_);
  else if (event->key() == Qt::Key_S)   // Backward
    scene.camera.MoveForward(-stepSize_);
  else if (event->key() == Qt::Key_D)   // Right
    scene.camera.MoveRight(stepSize_);
  else if (event->key() == Qt::Key_X)   // Up
    scene.camera.MoveUpward(stepSize_);
  else if (event->key() == Qt::Key_Z)   // Down
    scene.camera.MoveUpward(-stepSize_);
  else if (event->key() == Qt::Key_P)
    scene.LoadDefaultScene();
  else
  {
    changed = false;
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

    int dx = lastPos.x() - pos.x();
    int dy = lastPos.y() - pos.y();

    scene.camera.RotateY(0.1 * dx);         // Horizontal
    scene.camera.RotateX(0.1 * dy);         // Vertical

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

  for (const Object& obj : scene.objects)
  {
    // Initialize queue
    std::queue<OctreeNode*> q;
    q.push(&obj.octree->root);

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
      OctreeNode* node = q.front();
      q.pop();

      BoundingBox bb = node->bb;

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

	  if (node->children == nullptr)
		  continue;

      for (int i = 0; i < 8; ++i)
        q.push(node->children[i]);
    }
  }

  glDisableClientState(GL_VERTEX_ARRAY);
}

//--------------------------------------------------------------------------------
void GLWidget::drawCameraRay() const
{
  glBegin(GL_LINES);

  Triangle minTri;
  double minTime = std::numeric_limits<double>::max();

  Triangle tri;
  double time;
  bool hit = false;

  for (const Object& obj : scene.objects)
  {
    if (obj.octree->Query(debugRay_, tri, time) && time < minTime)
    {
      hit = true;
      minTri = tri;
      minTime = time;
    }
  }

  if (hit)
  {
    Vector3D point = debugRay_.Origin + minTime * debugRay_.Direction;
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

  for (const Object& obj : scene.objects)
  {
    for (const Triangle& triangle : obj.triangles)
    {
      for (const Vertex& vertex : triangle.Vertices)
      {
        glColor3f(obj.material.color.R, obj.material.color.G, obj.material.color.B);
        glNormal3f(vertex.Normal.X, vertex.Normal.Y, vertex.Normal.Z);
        glVertex3f(vertex.Position.X, vertex.Position.Y, vertex.Position.Z);
      }
    }
  }

  glEnd();
}