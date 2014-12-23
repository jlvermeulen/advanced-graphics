#define _USE_MATH_DEFINES

#include <glWidget.h>

#include <fstream>
#include <gl/GLU.h>
#include <math.h>
#include <ObjReader.h>
#include <QKeyEvent>
#include <QMouseEvent>
#include <Vector3D.h>
#include <queue>

//--------------------------------------------------------------------------------
GLWidget::GLWidget(QWidget* parent)
  : QGLWidget(parent),
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
  octree = Octree(triangles, 100, 5);
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

  // Draw triangles
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

  // Draw octree raster
  glDisable(GL_LIGHTING);
  glBegin(GL_LINES);

  glColor3f(1, 0, 0);

  std::queue<OctreeNode> q;
  q.push(octree.root);

  while (!q.empty())
  {
	  OctreeNode node = q.front();
	  q.pop();

	  BoundingBox bb = node.bb;
	  Vector3D verts[8];
	  for (int i = 0; i < 8; i++)
	  {
		  double x = i & 4 ? bb.Halfsize.X : -bb.Halfsize.X;
		  double y = i & 2 ? bb.Halfsize.Y : -bb.Halfsize.Y;
		  double z = i & 1 ? bb.Halfsize.Z : -bb.Halfsize.Z;

		  verts[i] = bb.Center + Vector3D(x, y, z);
	  }

	  drawLine(verts[0], verts[1]);
	  drawLine(verts[0], verts[2]);
	  drawLine(verts[0], verts[4]);
	  drawLine(verts[1], verts[3]);
	  drawLine(verts[1], verts[5]);
	  drawLine(verts[2], verts[3]);
	  drawLine(verts[2], verts[6]);
	  drawLine(verts[3], verts[7]);
	  drawLine(verts[4], verts[5]);
	  drawLine(verts[4], verts[6]);
	  drawLine(verts[5], verts[7]);
	  drawLine(verts[6], verts[7]);

	  for (OctreeNode& n : node.children)
		  q.push(n);
  }

  glEnd();
}

//--------------------------------------------------------------------------------
void GLWidget::keyPressEvent(QKeyEvent* event)
{
  float step = 0.1f;
  bool changed = false;

  if (event->key() == Qt::Key_W)        // Forward
  {
    changed = true;
    camera_.MoveForward(step);
  }
  else if (event->key() == Qt::Key_A)   // Left
  {
    changed = true;
    camera_.MoveRight(-step);
  }
  else if (event->key() == Qt::Key_S)   // Backward
  {
    changed = true;
    camera_.MoveForward(-step);
  }
  else if (event->key() == Qt::Key_D)   // Right
  {
    changed = true;
    camera_.MoveRight(step);
  }
  else if (event->key() == Qt::Key_X)   // Up
  {
    changed = true;
    camera_.MoveUpward(step);
  }
  else if (event->key() == Qt::Key_Z)   // Down
  {
    changed = true;
    camera_.MoveUpward(-step);
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
void GLWidget::drawLine(const Vector3D& v1, const Vector3D& v2) const
{
  glVertex3f(v1.X, v1.Y, v1.Z);
  glVertex3f(v2.X, v2.Y, v2.Z);
}